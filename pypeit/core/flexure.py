""" Module for flexure routines

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import copy
import inspect
import pathlib

import numpy as np

from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib

from astropy import stats
from astropy import units
from astropy.io import ascii
import scipy.signal
import scipy.optimize as opt
from scipy import interpolate

from linetools.spectra import xspectrum1d

from pypeit import msgs
from pypeit import utils
from pypeit.display import display
from pypeit.core.wavecal import autoid
from pypeit.core import arc
from pypeit.core import extract
from pypeit.core import fitting
from pypeit.core import qa
from pypeit.datamodel import DataContainer
from pypeit.images.detector_container import DetectorContainer
from pypeit.images.mosaic import Mosaic
from pypeit import specobj, specobjs
from pypeit import data

from IPython import embed


def spat_flexure_shift(sciimg, slits, debug=False, maxlag=20):
    """
    Calculate a rigid flexure shift in the spatial dimension
    between the slitmask and the science image.

    It is *important* to use original=True when defining the
    slitmask as everything should be relative to the initial slits

    Otherwise, the WaveTilts could get out of sync with science images

    Args:
        sciimg (`numpy.ndarray`_):
        slits (:class:`pypeit.slittrace.SlitTraceSet`):
        maxlag (:obj:`int`, optional):
            Maximum flexure searched for

    Returns:
        float:  The spatial flexure shift relative to the initial slits

    """
    # Mask -- Includes short slits and those excluded by the user (e.g. ['rdx']['slitspatnum'])
    slitmask = slits.slit_img(initial=True, exclude_flag=slits.bitmask.exclude_for_flexure)

    _sciimg = sciimg if slitmask.shape == sciimg.shape \
                else arc.resize_mask2arc(slitmask.shape, sciimg) 
    onslits = slitmask > -1
    corr_slits = onslits.astype(float).flatten()

    # Compute
    mean_sci, med_sci, stddev_sci = stats.sigma_clipped_stats(_sciimg[onslits])
    thresh =  med_sci + 5.0*stddev_sci
    corr_sci = np.fmin(_sciimg.flatten(), thresh)

    lags, xcorr = utils.cross_correlate(corr_sci, corr_slits, maxlag)
    xcorr_denom = np.sqrt(np.sum(corr_sci*corr_sci)*np.sum(corr_slits*corr_slits))
    xcorr_norm = xcorr / xcorr_denom
    # TODO -- Generate a QA plot
    tampl_true, tampl, pix_max, twid, centerr, ww, arc_cont, nsig \
            = arc.detect_lines(xcorr_norm, sigdetect=3.0, fit_frac_fwhm=1.5, fwhm=5.0,
                               cont_frac_fwhm=1.0, cont_samp=30, nfind=1, debug=debug)
    # No peak? -- e.g. data fills the entire detector
    if len(tampl) == 0:
        msgs.warn('No peak found in spatial flexure.  Assuming there is none..')
#        if debug:
#            embed(header='68 of flexure')
        return 0.

    # Find the peak
    xcorr_max = np.interp(pix_max, np.arange(lags.shape[0]), xcorr_norm)
    lag_max = np.interp(pix_max, np.arange(lags.shape[0]), lags)
    msgs.info('Spatial flexure measured: {}'.format(lag_max[0]))

    if debug:
        plt.figure(figsize=(14, 6))
        plt.plot(lags, xcorr_norm, color='black', drawstyle='steps-mid', lw=3, label='x-corr', linewidth=1.0)
        plt.plot(lag_max[0], xcorr_max[0], 'g+', markersize=6.0, label='peak')
        plt.title('Best shift = {:5.3f}'.format(lag_max[0]) + ',  corr_max = {:5.3f}'.format(xcorr_max[0]))
        plt.legend()
        plt.show()

    #tslits_shift = trace_slits.shift_slits(tslits_dict, lag_max)
    # Now translate the tilts

    #slitmask_shift = pixels.tslits2mask(tslits_shift)
    #slitmask_shift = slits.slit_img(flexure=lag_max[0])
    if debug:
        # Now translate the slits in the tslits_dict
        all_left_flexure, all_right_flexure, mask = slits.select_edges(flexure=lag_max[0])
        gpm = mask == 0
        viewer, ch = display.show_image(_sciimg)
        #display.show_slits(viewer, ch, left_flexure[:,gpm], right_flexure)[:,gpm]#, slits.id) #, args.det)
        #embed(header='83 of flexure.py')

    return lag_max[0]


def spec_flex_shift(obj_skyspec, arx_skyspec, arx_fwhm_pix, spec_fwhm_pix=None, mxshft=20, excess_shft="crash",
                    method="boxcar"):
    """ Calculate shift between object sky spectrum and archive sky spectrum

    Args:
        obj_skyspec (:class:`linetools.spectra.xspectrum1d.XSpectrum1d`):
            Spectrum of the sky related to our object
        arx_skyspec (:class:`linetools.spectra.xspectrum1d.XSpectrum1d`):
            Archived sky spectrum
        arx_fwhm_pix (:obj:`float`):
            Spectral FWHM (in pixels) of the archived sky spectrum.
        spec_fwhm_pix (:obj:`float`, optional):
            Spectral FWHM (in pixels) of the sky spectrum related to our object/slit.
        mxshft (:obj:`int`, optional):
            Maximum allowed shift from flexure;  note there are cases that
            have been known to exceed even 30 pixels.
        excess_shft (:obj:`str`, optional):
            Behavior of the code when a measured flexure exceeds ``mxshft``.
            Options are "crash", "set_to_zero", and "continue", where
            "set_to_zero" sets the shift to zero and moves on, and
            "continue" simply uses the large flexure shift value.
        method (:obj:`str`, optional):
            Which method is used for the spectral flexure correction.
            Two methods are available: 'boxcar' and 'slitcen' (see spec_flexure_slit()).
            In this routine, 'method' is only passed to final dict.

    Returns:
        dict: Contains flexure info.  Keys are:

          - polyfit= fit to the cross-correlation
          - shift= best shift in pixels
          - subpix= subpixelation of input spectrum
          - corr= correlation function
          - sky_spec= object sky spectrum used (rebinned, etc.)
          - arx_spec= archived sky spectrum used
          - corr_cen= center of the correlation function
          - smooth= Degree of smoothing of input spectrum to match archive
    """

    # TODO None of these routines should have dependencies on XSpectrum1d!

    msgs.warn("If we use Paranal, cut down on wavelength early on")

    # get gaussian sigma (pixels) for smoothing
    smooth_fwhm_pix = get_fwhm_gauss_smooth(arx_skyspec, obj_skyspec, arx_fwhm_pix, spec_fwhm_pix=spec_fwhm_pix)

    if smooth_fwhm_pix is None:
        # smooth_fwhm_pix is None if spec_fwhm_pix<0, i.e., the wavelength calibration is bad
        msgs.warn('No flexure correction could be computed for this slit/object')
        return None

    if smooth_fwhm_pix > 0:
        arx_skyspec = arx_skyspec.gauss_smooth(smooth_fwhm_pix)

    #Determine region of wavelength overlap
    min_wave = max(np.amin(arx_skyspec.wavelength.value), np.amin(obj_skyspec.wavelength.value))
    max_wave = min(np.amax(arx_skyspec.wavelength.value), np.amax(obj_skyspec.wavelength.value))

    #Smooth higher resolution spectrum by smooth_sig (flux is conserved!)
#    if np.median(obj_res) >= np.median(arx_res):
#        msgs.warn("New Sky has higher resolution than Archive.  Not smoothing")
        #obj_sky_newflux = ndimage.gaussian_filter(obj_sky.flux, smooth_sig)
#    else:
        #tmp = ndimage.gaussian_filter(arx_sky.flux, smooth_sig)
#        arx_skyspec = arx_skyspec.gauss_smooth(smooth_sig_pix*2*np.sqrt(2*np.log(2)))
        #arx_sky.flux = ndimage.gaussian_filter(arx_sky.flux, smooth_sig)

    # Define wavelengths of overlapping spectra
    keep_idx = np.where((obj_skyspec.wavelength.value>=min_wave) &
                         (obj_skyspec.wavelength.value<=max_wave))[0]
    #keep_wave = [i for i in obj_sky.wavelength.value if i>=min_wave if i<=max_wave]

    #Rebin both spectra onto overlapped wavelength range
    if len(keep_idx) <= 50:
        msgs.warn("Not enough overlap between sky spectra")
        return None

    # rebin onto object ALWAYS
    keep_wave = obj_skyspec.wavelength[keep_idx]
    arx_skyspec = arx_skyspec.rebin(keep_wave)
    obj_skyspec = obj_skyspec.rebin(keep_wave)
    # Trim edges (rebinning is junk there)
    arx_skyspec.data['flux'][0,:2] = 0.
    arx_skyspec.data['flux'][0,-2:] = 0.
    obj_skyspec.data['flux'][0,:2] = 0.
    obj_skyspec.data['flux'][0,-2:] = 0.

    # Set minimum to 0.  For bad rebinning and for pernicious extractions
    obj_skyspec.data['flux'][0,:] = np.maximum(obj_skyspec.data['flux'][0,:], 0.)
    arx_skyspec.data['flux'][0,:] = np.maximum(arx_skyspec.data['flux'][0,:], 0.)

    # Normalize spectra to unit average sky count
    norm = np.sum(obj_skyspec.flux.value)/obj_skyspec.npix
    norm2 = np.sum(arx_skyspec.flux.value)/arx_skyspec.npix
    if norm <= 0:
        msgs.warn("Bad normalization of object in flexure algorithm")
        msgs.warn("Will try the median")
        norm = np.median(obj_skyspec.flux.value)
        if norm <= 0:
            msgs.warn("Improper sky spectrum for flexure.  Is it too faint??")
            return None
    if norm2 <= 0:
        msgs.warn('Bad normalization of archive in flexure. You are probably using wavelengths '
                   'well beyond the archive.')
        return None
    obj_skyspec.flux = obj_skyspec.flux / norm
    arx_skyspec.flux = arx_skyspec.flux / norm2

    # Deal with bad pixels
    msgs.work("Need to mask bad pixels")

    # Deal with underlying continuum
    msgs.work("Consider taking median first [5 pixel]")
    everyn = obj_skyspec.npix // 20
    pypeitFit_obj, _ = fitting.iterfit(obj_skyspec.wavelength.value, obj_skyspec.flux.value,
                                       nord = 3,  kwargs_bspline={'everyn': everyn}, kwargs_reject={'groupbadpix':True,'maxrej':1},
                                       maxiter = 15, upper = 3.0, lower = 3.0)
    obj_sky_cont, _ = pypeitFit_obj.value(obj_skyspec.wavelength.value)

    obj_sky_flux = obj_skyspec.flux.value - obj_sky_cont
    pypeitFit_sky, _ = fitting.iterfit(arx_skyspec.wavelength.value, arx_skyspec.flux.value,
                                       nord = 3,  kwargs_bspline={'everyn': everyn}, kwargs_reject={'groupbadpix':True,'maxrej':1},
                                       maxiter = 15, upper = 3.0, lower = 3.0)
    arx_sky_cont, _ = pypeitFit_sky.value(arx_skyspec.wavelength.value)
    arx_sky_flux = arx_skyspec.flux.value - arx_sky_cont

    # Consider sharpness filtering (e.g. LowRedux)
    msgs.work("Consider taking median first [5 pixel]")

    #Cross correlation of spectra
    #corr = np.correlate(arx_skyspec.flux, obj_skyspec.flux, "same")
    corr = np.correlate(arx_sky_flux, obj_sky_flux, "same")

    #Create array around the max of the correlation function for fitting for subpixel max
    # Restrict to pixels within maxshift of zero lag
    lag0 = corr.size//2
    #mxshft = settings.argflag['reduce']['flexure']['maxshift']
    max_corr = np.argmax(corr[lag0-mxshft:lag0+mxshft]) + lag0-mxshft
    subpix_grid = np.linspace(max_corr-3., max_corr+3., 7)

    #Fit a 2-degree polynomial to peak of correlation function. JFH added this if/else to not crash for bad slits
    if np.any(np.isfinite(corr[subpix_grid.astype(int)])):
        fit = fitting.PypeItFit(xval=subpix_grid, yval=corr[subpix_grid.astype(int)],
                                func='polynomial', order=np.atleast_1d(2))
        fit.fit()
        max_fit = -0.5 * fit.fitc[1] / fit.fitc[2]

        shift = float(max_fit) - lag0
        # Deal with the case of shifts greater than ``mxshft``
        # We need to compare the absolute value of shift to ``mxshft``, since shift can be
        # positive or negative, while ``mxshft`` is generally only positive
        # We use the int of abs(shift) to avoid to trigger the error/warning for differences <1pixel
        # TODO :: I'm not convinced that we need int here...
        if int(abs(shift)) > mxshft:
            msgs.warn(f"Computed shift {shift:.1f} pix is "
                      f"larger than specified maximum {mxshft} pix.")

            if excess_shft == "crash":
                msgs.error(f"Flexure compensation failed for one of your{msgs.newline()}"
                           f"objects.  Either adjust the \"spec_maxshift\"{msgs.newline()}"
                           f"FlexurePar Keyword, or see the flexure documentation{msgs.newline()}"
                           f"for information on how to bypass this error using the{msgs.newline()}"
                           f"\"excessive_shift\" keyword.{msgs.newline()}"
                           "https://pypeit.readthedocs.io/en/release/flexure.html")

            elif excess_shft == "set_to_zero":
                msgs.warn("Flexure compensation failed for one of your objects.")
                msgs.warn("Setting the flexure correction shift to 0 pixels.")
                # Return the usual dictionary, but with a shift == 0
                shift = 0.0

            elif excess_shft == "continue":
                msgs.warn("Applying flexure shift larger than specified max!")

            elif excess_shft == "use_median":
                msgs.warn("Will try to use a flexure shift from other slit/object. "
                          "If not available, flexure correction will not be applied.")
                return None

            else:
                msgs.error(f"FlexurePar Keyword excessive_shift = \"{excess_shft}\" "
                           "not recognized.")

    else:
        fit = fitting.PypeItFit(xval=subpix_grid, yval=0.0*subpix_grid,
                                func='polynomial', order=np.atleast_1d(2))
        fit.fit()
        msgs.warn('Flexure compensation failed for one of your objects')
        return None

    #Calculate and apply shift in wavelength
    # shift = float(max_fit)-lag0
    msgs.info(f"Flexure correction of {shift:.3f} pixels")
    #model = (fit[2]*(subpix_grid**2.))+(fit[1]*subpix_grid)+fit[0]

    return dict(polyfit=fit, shift=shift, subpix=subpix_grid,
                corr=corr[subpix_grid.astype(int)], sky_spec=obj_skyspec, arx_spec=arx_skyspec,
                corr_cen=lag0, smooth=smooth_fwhm_pix, method=method)


def get_fwhm_gauss_smooth(arx_skyspec, obj_skyspec, arx_fwhm_pix, spec_fwhm_pix=None):
    """

    Args:
        arx_skyspec (:class:`linetools.spectra.xspectrum1d.XSpectrum1d`):
            Archived sky spectrum.
        obj_skyspec (:class:`linetools.spectra.xspectrum1d.XSpectrum1d`):
            Sky spectrum associated with the science target.
        arx_fwhm_pix (:obj:`float`):
            Spectral FWHM (in pixels) of the archived sky spectrum.
        spec_fwhm_pix (:obj:`float`, optional):
            Spectral FWHM (in pixels) of the sky spectrum related to our object.

    Returns:
        :obj:`float`: FWHM of the smoothing Gaussian in pixels.
    """
    # determine object spectral FWHM (in Angstrom) using obj_skyspec
    # if spec_fwhm_pix (typically from wave calibration) is None
    if spec_fwhm_pix is None:
        # pixels
        spec_fwhm_pix = autoid.measure_fwhm(obj_skyspec.flux.value, sigdetect=4., fwhm=4.)
        msgs.info('Measuring spectral FWHM using the boxcar extracted sky spectrum.')
        if spec_fwhm_pix is None:
            msgs.warn('Failed to measure the spectral FWHM using the boxcar extracted sky spectrum. '
                      'Not enough sky lines detected.')
            return None
    # object sky spectral dispersion (Angstrom/pixel)
    obj_disp = np.median(np.diff(obj_skyspec.wavelength.value))
    # Angstrom
    spec_fwhm = spec_fwhm_pix * obj_disp

    # determine arxiv sky spectral FWHM (in Angstrom)
    # arxiv sky spectral dispersion (Angstrom/pixel)
    arx_disp = np.median(np.diff(arx_skyspec.wavelength.value))
    arx_fwhm = arx_fwhm_pix * arx_disp

    msgs.info(f"Resolution (FWHM) of Archive={arx_fwhm:.2f} Ang and Observation={spec_fwhm:.2f} Ang")

    if spec_fwhm <= 0:
        msgs.warn('Negative spectral FWHM, likely due to a bad wavelength calibration.')
        return None

    # Determine fwhm of the smoothing gaussian
    # object sky spectral fwhm (Angstrom)
    obj_med_fwhm2 = np.power(spec_fwhm, 2)
    # arxiv sky spectral fwhm (Angstrom)
    arx_med_fwhm2 = np.power(arx_fwhm, 2)

    if obj_med_fwhm2 >= arx_med_fwhm2:
        smooth_fwhm = np.sqrt(obj_med_fwhm2-arx_med_fwhm2)  # Ang
        smooth_fwhm_pix = smooth_fwhm / arx_disp
    else:
        msgs.warn("Prefer archival sky spectrum to have higher resolution")
        smooth_fwhm_pix = 0.
        msgs.warn("New Sky has higher resolution than Archive.  Not smoothing")

    return smooth_fwhm_pix


def flexure_interp(shift, wave):
    """
    Perform interpolation on wave given a shift in pixels

    Args:
        shift (float):
            Shift in pixels
        wave (`numpy.ndarray`_):
            extracted wave of size nspec
        wavein (`numpy.ndarray`_, optional):
            Apply the shift to this array of wavelengths
    Returns:
        `numpy.ndarray`_: Wavelength scale corrected for spectral flexure

    """
    npix = wave.size
    x = np.linspace(0., 1., npix)
    f = interpolate.interp1d(x, wave, bounds_error=False, fill_value="extrapolate")
    twave = f(x + shift / (npix - 1))
    return twave


def spec_flex_shift_global(slit_specs, islit, sky_spectrum, arx_fwhm_pix, empty_flex_dict,
                           return_later_slits, flex_list, keys_to_update, spec_fwhm_pix=None,
                           mxshft=20, excess_shft="crash", method='slitcen'):
    """ Calculate flexure shifts using the sky spectrum extracted at the center of the slit

    Args:
        slit_specs (:obj:`list`):
            A list of `linetools.xspectrum1d`, one for each slit. The spectra stored in
            this list are sky spectra, extracted from the center of each slit.
        islit (:obj:`int`):
            Index of the slit where the sky spectrum related to our object is.
        sky_spectrum (:class:`linetools.spectra.xspectrum1d.XSpectrum1d`):
            Archived sky spectrum.
        arx_fwhm_pix (:obj:`float`):
            Spectral FWHM (in pixels) of the archived sky spectrum.
        empty_flex_dict (:obj:`dict`):
            Empty dictionary to be filled with flexure results.
        return_later_slits (:obj:`list`):
            List of slit indexes that failed the shift calcultion and we want to come back to
            to assign a value from a different slit.
        flex_list (:obj:`list`):
            A list of :obj:`dict` objects containing flexure results of each slit.
        keys_to_update (:obj:`list`):
            List of flexure dictionary keys that we need to update.
        spec_fwhm_pix (:obj:`float`, optional):
            Spectral FWHM (in pixels) of the sky spectrum related to our object.
        mxshft (:obj:`int`, optional):
            Maximum allowed shift from flexure. Passed to spec_flex_shift().
        excess_shft (:obj:`str`, optional):
            Behavior of the code when a measured flexure exceeds ``mxshft``.
            Passed to spec_flex_shift()
        method (:obj:`str`, optional):
            Which method is used for the spectral flexure correction.
            Two methods are available: 'boxcar' and 'slitcen' (see spec_flexure_slit()).
            Passed to spec_flex_shift().

    Returns:
        :obj:`list`: A list of :obj:`dict` objects containing flexure
        results of each slit. This is filled with a basically empty
        dict if the shift calculation failed for the relevant slit.

    """

    # Reset the flexure dictionary
    flex_dict = copy.deepcopy(empty_flex_dict)

    # Calculate the shift
    fdict = spec_flex_shift(slit_specs[islit], sky_spectrum, arx_fwhm_pix, mxshft=mxshft, excess_shft=excess_shft,
                            spec_fwhm_pix=spec_fwhm_pix, method=method)

    # Was it successful?
    if fdict is not None:
        # Update dict
        for key in keys_to_update[:-1]:
            flex_dict[key].append(fdict[key])
        # Interpolate
        sky_wave_new = flexure_interp(fdict['shift'], slit_specs[islit].wavelength.value)
        flex_dict['sky_spec'].append(xspectrum1d.XSpectrum1D.from_tuple((sky_wave_new, slit_specs[islit].flux.value)))
    else:
        # No success, come back to it later
        return_later_slits.append(islit)
        msgs.warn("Flexure shift calculation failed for this slit.")
        msgs.info("Will come back to this slit to attempt "
                  "to use saved estimates from other slits")

    # Append flex_dict, which will be an empty dictionary if the flexure failed for the all the slits
    flex_list.append(flex_dict.copy())
    return flex_list


def spec_flex_shift_local(slits, slitord, specobjs, islit, sky_spectrum, arx_fwhm_pix, empty_flex_dict,
                          return_later_slits, flex_list, keys_to_update, spec_fwhm_pix=None, mxshft=20,
                          excess_shft="crash", method='boxcar'):
    """ Calculate flexure shifts using the sky spectrum boxcar-extracted at the location of the detected objects

    Args:
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            Slit trace set.
        slitord (`numpy.ndarray`_):
            Array of slit/order numbers.
        specobjs (:class:`~pypeit.specobjs.Specobjs`, optional):
            Spectral extractions.
        islit (:obj:`int`):
            Index of the slit where the sky spectrum related to our object is.
        sky_spectrum (:class:`linetools.spectra.xspectrum1d.XSpectrum1d`):
            Archived sky spectrum.
        arx_fwhm_pix (:obj:`float`):
            Spectral FWHM (in pixels) of the archived sky spectrum.
        empty_flex_dict (:obj:`dict`):
            Empty dictionary to be filled with flexure results.
        return_later_slits (:obj:`list`):
            List of slit indexes that failed the shift calcultion and we want to come back to
            to assign a value from a different slit.
        flex_list (:obj:`list`):
            A list of :obj:`dict` objects containing flexure results of each slit.
        keys_to_update (:obj:`list`):
            List of flexure dictionary keys that we need to update.
        spec_fwhm_pix (:obj:`float`, optional):
            Spectral FWHM (in pixels) of the sky spectrum related to our object.
        mxshft (:obj:`int`, optional):
            Maximum allowed shift from flexure. Passed to spec_flex_shift().
        excess_shft (:obj:`str`, optional):
            Behavior of the code when a measured flexure exceeds ``mxshft``.
            Passed to spec_flex_shift()
        method (:obj:`str`, optional):
            Which method is used for the spectral flexure correction.
            Two methods are available: 'boxcar' and 'slitcen' (see spec_flexure_slit()).
            Passed to spec_flex_shift().
    Returns:
        :obj:`list`: A list of :obj:`dict` objects containing flexure
        results of each slit. This is filled with a basically empty
        dict if the shift calculation failed for the relevant slit.
    """

    # Reset the flexure dictionary
    flex_dict = copy.deepcopy(empty_flex_dict)

    # get objects in this slit
    i_slitord = slitord[islit]
    indx = specobjs.slitorder_indices(i_slitord)
    this_specobjs = specobjs[indx]

    # if no objects in this slit, append an empty dict
    if len(this_specobjs) == 0:
        msgs.info('No object extracted in this slit.')
        flex_list.append(empty_flex_dict.copy())
        return flex_list

    # Objects in this slit that failed and we want to come back to
    # to assign values from other objects in the same slit (if available)
    return_later_sobjs = []
    # Loop through objects
    for ss, sobj in enumerate(this_specobjs):
        if sobj is None or sobj['BOX_WAVE'] is None:  # Nothing extracted; only the trace exists
            msgs.info(f'Object # {ss} was not extracted.')
            # Update dict
            for key in keys_to_update:
                # append None
                flex_dict[key].append(None)
            continue
        msgs.info(f"Working on spectral flexure for object # {ss} in slit {slits.spat_id[islit]}")

        # get 1D spectrum for this object
        obj_sky = xspectrum1d.XSpectrum1D.from_tuple((sobj.BOX_WAVE, sobj.BOX_COUNTS_SKY))

        # Calculate the shift
        fdict = spec_flex_shift(obj_sky, sky_spectrum, arx_fwhm_pix, mxshft=mxshft, excess_shft=excess_shft,
                                spec_fwhm_pix=spec_fwhm_pix, method=method)

        if fdict is not None:
            # Update dict
            for key in keys_to_update:
                flex_dict[key].append(fdict[key])
        else:
            # No success, come back to it later
            return_later_sobjs.append(ss)
            msgs.warn("Flexure shift calculation failed for this spectrum.")
            msgs.info("Will come back to this spectrum to attempt "
                      "to use saved estimates from other slits/objects")

    # Check if we need to go back
    if (len(return_later_sobjs) > 0) and (len(flex_dict['shift']) > 0):
        msgs.warn(f'Flexure shift calculation failed for {len(return_later_sobjs)} '
                  f'object(s) in slit {slits.spat_id[islit]}')
        # get the median shift among all objects in this slit
        idx_med_shift = np.where(flex_dict['shift'] == np.percentile(flex_dict['shift'], 50,
                                                                     interpolation='nearest'))[0][0]
        msgs.info(f"Median value of the measured flexure shifts in this slit, equal to "
                  f"{flex_dict['shift'][idx_med_shift]:.3f} pixels, will be used")

        # assign the median shift to the failed objects
        for obj_idx in return_later_sobjs:
            # Update dict
            for key in keys_to_update[:-1]:
                # insert the median value at the location of the object that failed the calculation
                flex_dict[key].insert(obj_idx, flex_dict[key][idx_med_shift])
            # Interpolate
            sky_wave_new = flexure_interp(flex_dict['shift'][obj_idx], this_specobjs[obj_idx].BOX_WAVE)
            flex_dict['sky_spec'].insert(obj_idx, xspectrum1d.XSpectrum1D.from_tuple(
                (sky_wave_new, this_specobjs[obj_idx].BOX_COUNTS_SKY)))

    # if flexure failed for every objects in this slit, save for later to use value from other slits
    elif (len(return_later_sobjs) > 0) and (len(flex_dict['shift']) == 0):
        return_later_slits.append(islit)

    # Append flex_dict, which will be an empty dictionary if the flexure failed for the whole slit
    flex_list.append(flex_dict.copy())

    return flex_list


def spec_flexure_slit(slits, slitord, slit_bpm, sky_file, method="boxcar", specobjs=None,
                      slit_specs=None, wv_calib=None, mxshft=None, excess_shft="crash"):
    """Calculate the spectral flexure for every slit (global) or object (local)

    Args:
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            Slit trace set
        slitord (`numpy.ndarray`_):
            Array of slit/order numbers
        slit_bpm (`numpy.ndarray`_):
            True = masked slit
        sky_file (:obj:`str`):
            Sky file
        method (:obj:`str`, optional):
            Two methods are available:
                - 'boxcar': Recommended for object extractions. This
                  method uses the boxcar extracted sky and wavelength
                  spectra from the input specobjs
                - 'slitcen': Recommended when no objects are being
                  extracted. This method uses a spectrum (stored in
                  slitspecs) that is extracted from the center of
                  each slit.
        specobjs (:class:`~pypeit.specobjs.Specobjs`, optional):
            Spectral extractions
        slit_specs (:obj:`list`, optional):
            A list of linetools.xspectrum1d, one for each slit. The spectra stored in
            this list are sky spectra, extracted from the center of each slit.
        wv_calib (:class:`pypeit.wavecalib.WaveCalib`):
            Wavelength calibration object
        mxshft (int, optional):
            Passed to spec_flex_shift()
        excess_shft (str, optional):
            Passed to spec_flex_shift()

    Returns:
        :obj:`list`: A list of :obj:`dict` objects containing flexure
        results of each slit. This is filled with a basically empty
        dict if the slit is skipped.
    """
    msgs.work("Consider doing 2 passes in flexure as in LowRedux")

    # Determine the method
    slit_cen = True if (specobjs is None) or (method == "slitcen") else False

    # Load Archival sky spectrum
    sky_spectrum, arx_fwhm_pix = get_archive_spectrum(sky_file)

    # Initialise the flexure list for each slit
    flex_list = []

    # initiate list of slits to come back to if flexure calculation failed
    return_later_slits = []

    # empty dict
    empty_flex_dict = dict(polyfit=[], shift=[], subpix=[], corr=[],
                           corr_cen=[], spec_file=sky_file, smooth=[],
                           arx_spec=[], sky_spec=[], method=[])

    # flex dict keys that we need to update through the routine
    keys_to_update = ['polyfit', 'shift', 'subpix', 'corr', 'corr_cen', 'smooth', 'method', 'arx_spec', 'sky_spec']

    # Loop over slits
    # good slits
    gdslits = np.where(np.logical_not(slit_bpm))[0]

    for islit in range(slits.nslits):
        msgs.info(f"Working on spectral flexure of slit: {slits.spat_id[islit]}")

        # If no objects on this slit append an empty dictionary
        if islit not in gdslits:
            flex_list.append(empty_flex_dict.copy())
            continue

        # get spectral FWHM (in pixels) if available
        spec_fwhm_pix = None
        if wv_calib is not None:
            iwv = np.where(wv_calib.spat_ids == slits.spat_id[islit])[0][0]
            # Allow for wavelength failures
            if wv_calib.wv_fits is not None and wv_calib.wv_fits[iwv].fwhm is not None:
                spec_fwhm_pix = wv_calib.wv_fits[iwv].fwhm

        if slit_cen:
            # global flexure
            flex_list = spec_flex_shift_global(slit_specs, islit, sky_spectrum, arx_fwhm_pix, empty_flex_dict,
                                               return_later_slits, flex_list, keys_to_update,
                                               spec_fwhm_pix=spec_fwhm_pix, mxshft=mxshft, excess_shft=excess_shft)
        else:
            # local flexure
            flex_list = spec_flex_shift_local(slits, slitord, specobjs, islit, sky_spectrum, arx_fwhm_pix,
                                              empty_flex_dict, return_later_slits, flex_list, keys_to_update,
                                              spec_fwhm_pix=spec_fwhm_pix, mxshft=mxshft, excess_shft=excess_shft)

    # Check if we need to go back to some failed slits
    if len(return_later_slits) > 0:
        msgs.warn(f'Flexure shift calculation failed for {len(return_later_slits)} slits')
        # take the median value to deal with the cases when there are more than one shift per slit (e.g., local flexure)
        saved_shifts = np.array([np.percentile(flex['shift'], 50, interpolation='nearest')
                                 if len(flex['shift']) > 0 else None for flex in flex_list])
        if np.all(saved_shifts == None):
            # If all the elements in saved_shifts are None means that there are no saved shifts available
            msgs.warn(f'No previously saved flexure shift estimates available. '
                      f'Flexure corrections cannot be performed.')
            for islit in range(slits.nslits):
                # we append an empty dictionary
                flex_list.append(empty_flex_dict.copy())
        else:
            # get the median shift value among all slit
            med_shift = np.percentile(saved_shifts[saved_shifts!= None], 50, interpolation='nearest')
            # in which slit the median is?
            islit_med_shift = np.where(saved_shifts == med_shift)[0][0]
            msgs.info(f"Median value of all the measured flexure shifts, equal to "
                      f"{saved_shifts[islit_med_shift]:.3f} pixels, will be used")

            # global flexure
            if slit_cen:
                # get the dict where the med shift is
                fdict = flex_list[islit_med_shift].copy()

                # assign fdict to the failed slits
                for sidx in return_later_slits:
                    # Reset the dict
                    flex_dict = copy.deepcopy(empty_flex_dict)
                    # Update dict
                    for key in keys_to_update[:-1]:
                        flex_dict[key].append(fdict[key][0])
                    # Interpolate
                    sky_wave_new = flexure_interp(fdict['shift'][0], slit_specs[sidx].wavelength.value)
                    flex_dict['sky_spec'].append(xspectrum1d.XSpectrum1D.from_tuple((sky_wave_new, slit_specs[sidx].flux.value)))

                    # insert flex_dict in flex_list at the location of the slit that failed the calculation
                    flex_list[sidx] = flex_dict
            # local flexure
            else:
                # get the dict where the med shift is
                idx_med_shift = np.where(flex_list[islit_med_shift]['shift'] == med_shift)[0][0]
                fdict = copy.deepcopy(empty_flex_dict)
                for key in keys_to_update[:-1]:
                    fdict[key].append(flex_list[islit_med_shift][key][idx_med_shift])

                # assign fdict to the failed object
                for sidx in return_later_slits:
                    # get objects in this slit
                    i_slitord = slitord[sidx]
                    indx = specobjs.slitorder_indices(i_slitord)

                    for i in range(len(specobjs[indx])):
                        # Reset the dict
                        flex_dict = copy.deepcopy(empty_flex_dict)
                        # Update dict
                        for key in keys_to_update[:-1]:
                            flex_dict[key].append(fdict[key][0])
                        # Interpolate
                        sky_wave_new = flexure_interp(fdict['shift'][0], specobjs[indx][i].BOX_WAVE)
                        flex_dict['sky_spec'].append(
                            xspectrum1d.XSpectrum1D.from_tuple((sky_wave_new, specobjs[indx][i].BOX_COUNTS_SKY)))
                        # insert flex_dict in flex_list at the location of the slit that failed the calculation
                        flex_list[sidx] = flex_dict
    return flex_list


def spec_flexure_slit_global(sciImg, waveimg, global_sky, par, slits, slitmask, trace_spat, gd_slits, wv_calib, pypeline, det):
    """Calculate the spectral flexure for every slit

    Args:
        sciImg  (:class:`~pypeit.images.pypeitimage.PypeItImage`):
            Science image.
        waveimg (`numpy.ndarray`_):
            Wavelength image - shape (nspec, nspat)
        global_sky (`numpy.ndarray`_):
            2D array of the global_sky fit - shape (nspec, nspat)
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            Parameters of the reduction.
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            Slit trace set
        slitmask (`numpy.ndarray`_):
            An image with the slit index identified for each pixel (returned from slittrace.slit_img).
        trace_spat (`numpy.ndarray`_):
            Spatial pixel values (usually the center of each slit) where the sky spectrum will be extracted.
            The shape of this array should be (nspec, nslits)
        gd_slits (`numpy.ndarray`_):
            True = good slit
        wv_calib (:class:`pypeit.wavecalib.WaveCalib`):
            Wavelength calibration
        pypeline (:obj:`str`):
            Name of the ``PypeIt`` pipeline method.  Allowed options are
            MultiSlit, Echelle, or IFU.
        det (:obj:`str`):
            The name of the detector or mosaic from which the spectrum will be
            extracted.  For example, DET01.

    Returns:
        :obj:`list`: A list of :obj:`dict` objects containing flexure
        results of each slit. This is filled with a basically empty
        dict if the slit is skipped.
    """
    # TODO :: Need to think about spatial flexure - is the appropriate spatial flexure already included in trace_spat via left/right slits?
    slit_specs = []
    # get boxcar radius
    box_radius = par['reduce']['extraction']['boxcar_radius']
    for ss in range(slits.nslits):
        if not gd_slits[ss]:
            slit_specs.append(None)
            continue
        thismask = (slitmask == slits.spat_id[ss])
        inmask = sciImg.select_flag(invert=True) & thismask
        # Pack
        slit_specs.append(get_sky_spectrum(sciImg.image, sciImg.ivar, waveimg, inmask, global_sky,
                                           box_radius, slits, trace_spat[:, ss], pypeline, det))

    # Measure flexure
    flex_list = spec_flexure_slit(slits, slits.slitord_id, np.logical_not(gd_slits),
                                  par['flexure']['spectrum'],
                                  method=par['flexure']['spec_method'],
                                  mxshft=par['flexure']['spec_maxshift'],
                                  excess_shft=par['flexure']['excessive_shift'],
                                  specobjs=None, slit_specs=slit_specs, wv_calib=wv_calib)
    return flex_list


def get_archive_spectrum(sky_file):
    """ Load an archival sky spectrum

    Args:
        sky_file (:obj:`str`):
            Sky file

    Returns:
        (:obj:`XSpectrum1D`): Sky spectrum
        (float): FWHM of the sky lines in pixels.
    """
    # Load Archive. Save the fwhm to avoid the performance hit from calling it on the archive sky spectrum
    # multiple times
    sky_spectrum = data.load_sky_spectrum(sky_file)
    # get arxiv sky spectrum resolution (FWHM in pixels)
    arx_fwhm_pix = autoid.measure_fwhm(sky_spectrum.flux.value, sigdetect=4., fwhm=4.)
    if arx_fwhm_pix is None:
        msgs.error('Failed to measure the spectral FWHM of the archived sky spectrum. '
                   'Not enough sky lines detected.')
    return sky_spectrum, arx_fwhm_pix


def get_sky_spectrum(sciimg, ivar, waveimg, thismask, global_sky, box_radius, slits, trace_spat, pypeline, det):
    """ Obtain a boxcar extraction of the sky spectrum

    Args:
        sciimg  (`numpy.ndarray`_):
            Science image - shape (nspec, nspat)
        ivar  (`numpy.ndarray`_):
            Inverse variance of the science image - shape (nspec, nspat)
        waveimg (`numpy.ndarray`_):
            Wavelength image - shape (nspec, nspat)
        thismask (`numpy.ndarray`_):
            Good pixel mask (True=good) that indicates the pixels that should be included in the boxcar extraction
        global_sky (`numpy.ndarray`_):
            2D array of the global_sky fit - shape (nspec, nspat)
        box_radius (float):
            Radius of the boxcar extraction (in pixels)
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            Slit trace set
        trace_spat (`numpy.ndarray`_):
            Spatial pixel values (usually the center of each slit) where the sky spectrum will be extracted.
            The shape of this array should be (nspec, nslits)
        pypeline (:obj:`str`):
            Name of the ``PypeIt`` pipeline method.  Allowed options are
            MultiSlit, Echelle, or IFU.
        det (:obj:`str`):
            The name of the detector or mosaic from which the spectrum will be
            extracted.  For example, DET01.

    Returns:
        (:obj:`XSpectrum1D`): Sky spectrum
    """
    spec = specobj.SpecObj(PYPELINE=pypeline, SLITID=-1, DET=str(det))
    spec.trace_spec = np.arange(slits.nspec)
    spec.TRACE_SPAT = trace_spat
    spec.BOX_RADIUS = box_radius
    # Extract
    extract.extract_boxcar(sciimg, ivar, thismask, waveimg, global_sky, spec)
    slit_wave, slit_sky = spec.BOX_WAVE, spec.BOX_COUNTS_SKY
    # TODO :: Need to remove this XSpectrum1D dependency - it is required in:  flexure.spec_flex_shift
    obj_skyspec = xspectrum1d.XSpectrum1D.from_tuple((slit_wave, slit_sky))
    return obj_skyspec


def spec_flexure_corrQA(ax, this_flex_dict, cntr, name):
    # Fit
    fit = this_flex_dict['polyfit'][cntr]
    if fit is not None:
        xval = np.linspace(-10., 10, 100) + this_flex_dict['corr_cen'][cntr]  # + flex_dict['shift'][o]
        # model = (fit[2]*(xval**2.))+(fit[1]*xval)+fit[0]
        model = fit.eval(xval)
        # model = utils.func_val(fit, xval, 'polynomial')
        mxmod = np.max(model)
        ylim_min = np.min(model / mxmod) if np.isfinite(np.min(model / mxmod)) else 0.0
        ylim = [ylim_min, 1.3]
        ax.plot(xval - this_flex_dict['corr_cen'][cntr], model / mxmod, 'k-')
        # Measurements
        ax.scatter(this_flex_dict['subpix'][cntr] - this_flex_dict['corr_cen'][cntr],
                   this_flex_dict['corr'][cntr] / mxmod, marker='o')
        # Final shift
        ax.plot([this_flex_dict['shift'][cntr]] * 2, ylim, 'g:')
        # Label
        ax.text(0.5, 0.25, name, transform=ax.transAxes, size='large', ha='center')
        ax.text(0.5, 0.15, 'flex_shift = {:g}'.format(this_flex_dict['shift'][cntr]),
                transform=ax.transAxes, size='large', ha='center')  # , bbox={'facecolor':'white'})
        # Axes
        ax.set_ylim(ylim)
        ax.set_xlabel('Lag')
    else:
        ax.text(0.5, 0.25, name, transform=ax.transAxes, size='large', ha='center')
        ax.text(0.5, 0.15, 'flex_shift calculation failed', transform=ax.transAxes, size='large', ha='center')
        # Axes
        ax.set_xlabel('Lag')


def spec_flexure_qa(slitords, bpm, basename, flex_list,
                    specobjs=None, out_dir=None):
    """
    Generate QA for the spectral flexure calculation

    Args:
        slitords (`numpy.ndarray`_):
            Array of slit/order numbers
        bpm (`numpy.ndarray`_):
            Boolean mask; True = masked slit
        basename (str):
            Used to generate the output file name
        flex_list (list):
            list of :obj:`dict` objects containing the flexure information
        specobjs (:class:`~pypeit.specobjs.Specobjs`, optional):
            Spectrally extracted objects
        out_dir (str, optonal):
            Path to the output directory for the QA plots.  If None, the current
            is used.
    """
    plt.rcdefaults()
    plt.rcParams['font.family'] = 'serif'

    # What type of QA are we doing
    slit_cen = False
    if specobjs is None: slit_cen = True

    # Grab the named of the method
    method = inspect.stack()[0][3]

    # Mask
    gdslits = np.where(np.invert(bpm))[0]

    # Loop over slits, and then over objects here
    for islit in gdslits:
        # Slit/order number
        slitord = slitords[islit]

        this_flex_dict = flex_list[islit]
        # Check that the default was overwritten
        if len(this_flex_dict['shift']) == 0 or \
                (len(this_flex_dict['shift']) > 0 and np.all([ss is None for ss in this_flex_dict['shift']])):
            continue

        # Parse and Setup
        if slit_cen:
            nobj = 1
            ncol = 1
        else:
            indx = specobjs.slitorder_indices(slitord)
            this_specobjs = specobjs[indx]
            nobj = np.sum(indx)
            if nobj == 0:
                continue
            ncol = min(3, nobj)

        nrow = nobj // ncol + ((nobj % ncol) > 0)
        # Outfile, one QA file per slit
        outfile = qa.set_qa_filename(basename, method + '_corr', slit=slitord, out_dir=out_dir)
        plt.figure(figsize=(8, 5.0))
        plt.clf()
        gs = gridspec.GridSpec(nrow, ncol)
        # Correlation QA
        if slit_cen:
            ax = plt.subplot(gs[0, 0])
            spec_flexure_corrQA(ax, this_flex_dict, 0, 'Slit Center')
        else:
            iplt = 0
            for ss, specobj in enumerate(this_specobjs):
                if specobj is None or (specobj.BOX_WAVE is None and specobj.OPT_WAVE is None):
                    continue
                ax = plt.subplot(gs[iplt//ncol, iplt % ncol])
                spec_flexure_corrQA(ax, this_flex_dict, ss, '{:s}'.format(specobj.NAME))
                iplt += 1
        # Finish
        plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
        plt.savefig(outfile)#, dpi=400)
        plt.close()

        # Sky line QA (just one object)
        if slit_cen:
            iobj = 0
        else:
            # only show the first object in this slit that does not have None shift
            iobj = np.where([ss is not None for ss in this_flex_dict['shift']])[0][0]
            specobj = this_specobjs[iobj]

        # Repackage
        sky_spec = this_flex_dict['sky_spec'][iobj]
        arx_spec = this_flex_dict['arx_spec'][iobj]
        min_wave = max(np.amin(arx_spec.wavelength.value), np.amin(sky_spec.wavelength.value))*units.AA
        max_wave = min(np.amax(arx_spec.wavelength.value), np.amax(sky_spec.wavelength.value))*units.AA

        # Sky lines
        sky_lines = np.array([3370.0, 3914.0, 4046.56, 4358.34, 5577.338, 6300.304,
                              7340.885, 7993.332, 8430.174, 8919.610, 9439.660,
                              10013.99, 10372.88])*units.AA
        dwv = 20.*units.AA
        gdsky = np.where((sky_lines > min_wave) & (sky_lines < max_wave))[0]
        if len(gdsky) == 0:
            msgs.warn("No sky lines for Flexure QA")
            continue
        if len(gdsky) > 6:
            idx = np.array([0, 1, len(gdsky)//2, len(gdsky)//2+1, -2, -1])
            gdsky = gdsky[idx]

        # Outfile
        outfile = qa.set_qa_filename(basename, method+'_sky', slit=slitord, out_dir=out_dir)
        # Figure
        plt.figure(figsize=(8, 5.0))
        plt.clf()
        nrow, ncol = 2, 3
        gs = gridspec.GridSpec(nrow, ncol)
        if slit_cen:
            plt.suptitle('Sky Comparison for Slit Center', y=0.99)
        else:
            plt.suptitle('Sky Comparison for {:s}'.format(specobj.NAME), y=0.99)

        for ii, igdsky in enumerate(gdsky):
            skyline = sky_lines[igdsky]
            ax = plt.subplot(gs[ii//ncol, ii % ncol])
            # Norm
            pix1 = np.where(np.abs(sky_spec.wavelength-skyline) < dwv)[0]
            pix2 = np.where(np.abs(arx_spec.wavelength-skyline) < dwv)[0]
            f1 = np.sum(sky_spec.flux[pix1])
            f2 = np.sum(arx_spec.flux[pix2])
            norm = f1/f2
            # Plot
            ax.plot(sky_spec.wavelength[pix1], sky_spec.flux[pix1], 'k-', label='Obj',
                    drawstyle='steps-mid')
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
        plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
        plt.savefig(outfile)#, dpi=400)
        plt.close()
        msgs.info("Wrote spectral flexure QA: {}".format(outfile))

    plt.rcdefaults()


def calculate_image_phase(imref, imshift, gpm_ref=None, gpm_shift=None, maskval=None):
    """
    Perform a masked cross-correlation and optical flow calculation to robustly
    estimate the subpixel shifts of two images.

    If gpm_ref, gpm_shift, and maskval are all None, no pixels will be masked

    This routine (optionally) requires skimage to calculate the image phase. If
    skimage is not installed, a standard (unmasked) cross-correlation is used.


    Args:
        im_ref (`numpy.ndarray`_):
            Reference image
        imshift (`numpy.ndarray`_):
            Image that we want to measure the shift of (relative to im_ref)
        gpm_ref (`numpy.ndarray`_):
            Mask of good pixels (True = good) in the reference image
        gpm_shift (`numpy.ndarray`_):
            Mask of good pixels (True = good) in the shifted image
        maskval (float, optional):
            If gpm_ref and gpm_shift are both None, a single value can be specified
            and this value will be masked in both images.

    Returns:
        ra_diff (float):
            Relative shift (in pixels) of image relative to im_ref (x direction).
            In order to align image with im_ref, ra_diff should be added to the
            x-coordinates of image
        dec_diff (float):
            Relative shift (in pixels) of image relative to im_ref (y direction).
            In order to align image with im_ref, dec_diff should be added to the
            y-coordinates of image

    """
    # Do some checks first
    try:
        from skimage.registration import optical_flow_tvl1, phase_cross_correlation
    except ImportError:
        msgs.warn("scikit-image is not installed. Adopting a basic image cross-correlation")
        return calculate_image_offset(imref, imshift)
    if imref.shape != imshift.shape:
        msgs.warn("Input images shapes are not equal. Adopting a basic image cross-correlation")
        return calculate_image_offset(imref, imshift)
    # Set the masks
    if gpm_ref is None:
        gpm_ref = np.ones(imref.shape, dtype=bool) if maskval is None else imref != maskval
    if gpm_shift is None:
        gpm_shift = np.ones(imshift.shape, dtype=bool) if maskval is None else imshift != maskval
    # Get a crude estimate of the shift
    shift = phase_cross_correlation(imref, imshift, reference_mask=gpm_ref, moving_mask=gpm_shift).astype(int)
    # Extract the overlapping portion of the images
    exref = imref.copy()
    exshf = imshift.copy()
    if shift[0] != 0:
        if shift[0] < 0:
            exref = exref[:shift[0], :]
            exshf = exshf[-shift[0]:, :]
        else:
            exref = exref[shift[0]:, :]
            exshf = exshf[:-shift[0], :]
    if shift[1] != 0:
        if shift[1] < 0:
            exref = exref[:, :shift[1]]
            exshf = exshf[:, -shift[1]:]
        else:
            exref = exref[:, shift[1]:]
            exshf = exshf[:, :-shift[1]]
    # Compute the flow vector for a fine correction to the cross-correlation
    v, u = optical_flow_tvl1(exref, exshf)
    shift = shift.astype(float)
    shift[0] -= np.median(v)
    shift[1] -= np.median(u)
    # Return the total estimated shift
    return shift[0], shift[1]


def calculate_image_offset(im_ref, image, nfit=3):
    """Calculate the x,y offset between two images

    Args:
        im_ref (`numpy.ndarray`_):
            Reference image
        image (`numpy.ndarray`_):
            Image that we want to measure the shift of (relative to im_ref)
        nfit (int, optional):
            Number of pixels (left and right of the maximum) to include in
            fitting the peak of the cross correlation.

    Returns:
        tuple: Returns two floats, the x and y offset of the image.
          - ra_diff --  Relative shift (in pixels) of image relative to im_ref (x direction).
            In order to align image with im_ref, ra_diff should be added to the
            x-coordinates of image
          - dec_diff  -- Relative shift (in pixels) of image relative to im_ref (y direction).
            In order to align image with im_ref, dec_diff should be added to the
            y-coordinates of image
    """
    # Subtract median (should be close to zero, anyway)
    image -= np.median(image)
    im_ref -= np.median(im_ref)

    # cross correlate (note, convolving seems faster)
    ccorr = scipy.signal.correlate2d(im_ref, image, boundary='fill', mode='same')
    #ccorr = scipy.signal.fftconvolve(im_ref, image[::-1, ::-1], mode='same')

    # Find the maximum
    amax = np.unravel_index(np.argmax(ccorr), ccorr.shape)

    # Extract a small region around the maximum, and check the limits
    xlo, xhi = amax[0]-nfit, amax[0] + nfit+1
    ylo, yhi = amax[1]-nfit, amax[1] + nfit+1
    if xlo < 0: xlo = 0
    if xhi > ccorr.shape[0]-1: xhi = ccorr.shape[0]-1
    if ylo < 0: ylo = 0
    if yhi > ccorr.shape[1]-1: yhi = ccorr.shape[1]-1
    x = np.arange(xlo, xhi)
    y = np.arange(ylo, yhi)
    # Setup some initial parameters
    initial_guess = (np.max(ccorr), amax[0], amax[1], 3, 3, 0, 0)
    xx, yy = np.meshgrid(x, y, indexing='ij')

    # Fit the neighborhood of the maximum with a Gaussian to calculate the offset
    popt, _ = opt.curve_fit(fitting.twoD_Gaussian, (xx, yy), ccorr[xlo:xhi, ylo:yhi].ravel(), p0=initial_guess)
    # Return the RA and DEC shift, in pixels
    xoff = 1 - (ccorr.shape[0] % 2)  # Need to add 1 for even shaped array
    yoff = 1 - (ccorr.shape[1] % 2)  # Need to add 1 for even shaped array
    return xoff + popt[1] - ccorr.shape[0]//2, yoff+popt[2] - ccorr.shape[1]//2


def sky_em_residuals(wave:np.ndarray, flux:np.ndarray,
                     ivar:np.ndarray, sky_waves:np.ndarray,
                     plot=False, noff=5., nfit_min=20):
    """Calculate residuals and other metrics for a set of
    input sky emission lines 

    Args:
        wave (`numpy.ndarray`_):
            Wavelengths (in air!)
        flux (`numpy.ndarray`_):
            Fluxes
        ivar (`numpy.ndarray`_):
            Inverse variance
        sky_waves (`numpy.ndarray`_):
            Skyline wavelengths (in air!)
        plot (bool, optional):
            If true, plot the residuals
        noff (int, optional):
            Range in Ang to analyze labout emission line. Defaults to 5.
        nfit_min (int, optional):
            Minimum number of pixels required to do a fit. Defaults to 20.

    Returns:
        tuple of `numpy.ndarray`_ -- sky line wavelength of good lines, wavelength offset,
            error in wavelength offset, sky line width,
            error in sky line width
    """


    dwave = []
    diff = []
    diff_err  = []
    los = []
    los_err= []

    good_ivar = ivar > 0

    # Loop on known sky lines
    for line in sky_waves: 
        wline = [line-noff,line+noff] 
        mw    = (wave > wline[0]) & (wave < wline[1]) & good_ivar
        
        # Reuire minimum number
        if np.sum(mw) <= nfit_min:
            continue

        p=[0,0,0,0]
        # Guess
        p0 = list(fitting.guess_gauss(wave[mw], flux[mw]))
        # Fit
        try:
            p, pcov = fitting.fit_gauss(wave[mw], flux[mw], w_out=1./np.sqrt(ivar[mw]),
                                        guesses=p0, nparam=4)
        except RuntimeError as e:
            msgs.warn('First attempt at Gaussian fit failed, ending with RuntimeError.  Original '
                      f'exception: {e.args[0]}  Assuming this is because it hit the maximum '
                      'number of function evaluations.  Trying again with a maximum of 10000.')
            # Try again with larger limit on the number of function evaluations
            p, pcov = fitting.fit_gauss(wave[mw], flux[mw], w_out=1./np.sqrt(ivar[mw]),
                                        guesses=p0, nparam=4, maxfev=10000)

        perr = np.sqrt(np.diag(pcov))
        #except:
        #    p=p0
        #    p[2] = -99
        #    perr=p0

        # Continue
        d = p[2] - line

        # For debugging
        if plot:
            gfit = fitting.gauss_4deg(wave[mw],*p)
            plt.figure(figsize=(8,3)) 
            plt.plot(wave[mw],gfit,'g')
            plt.plot(wave[mw],flux[mw])
            plt.title('{} {:0.2f} diff= {:0.3f}'.format(line,p[3],d))
            plt.show()

        # Check
        if not np.isfinite(perr[2]):
            perr[2] = 1000.
        # Save
        dwave = np.append(dwave,line)
        diff = np.append(diff,d)
        diff_err = np.append(diff_err,perr[2])
        los = np.append(los,p[3])
        los_err = np.append(los_err,perr[3])

    # Cut on quality
    m=(diff_err < 0.1) & (diff_err > 0.0)
    # Return
    return dwave[m], diff[m], diff_err[m], los[m], los_err[m]


# TODO -- Consider separating the methods from the DataContainer as per calibrations
class MultiSlitFlexure(DataContainer):
    """
    Class to perform multi-detector flexure analysis.

    Based on code written by Marla Geha for DEIMOS.

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_multislitflexure.rst
    """

    # Set the version of this class
    version = '1.1.0'

    datamodel = {'s1dfile': dict(otype=str, descr='spec1d filename'), 
                 'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'ndet': dict(otype=int, descr='Number of detectors per spectrum'),
                 'nslits': dict(otype=int, descr='Number of slits'),
                 'is_msc': dict(otype=np.ndarray, atype=(int, np.integer),
                                descr='Flag that the "det" is the mosaic ID (ndet, nslits)'),
                 'det': dict(otype=np.ndarray, atype=(int, np.integer),
                             descr='Integer identifiers for the detector or mosaic (ndet, nslits)'),
                 'SN': dict(otype=np.ndarray, atype=np.floating, descr='S/N (ndet, nslits)'),
                 'slitid': dict(otype=np.ndarray, atype=np.floating, descr='Slit ID (nslits)'),
                 'mn_wv': dict(otype=np.ndarray, atype=np.floating,
                               descr='Mininum wavelength of the slit [Ang] (nslits)'),
                 'indiv_fit_slope': dict(otype=np.ndarray, atype=np.floating,
                                         descr='Fits to each slit individually (nslits)'),
                 'indiv_fit_b': dict(otype=np.ndarray, atype=np.floating,
                                     descr='Same as above but for b (nslits)'),
                 'indiv_fit_los': dict(otype=np.ndarray, atype=np.floating,
                                       descr='Same as above but for line width (nslits)'),
                 'fit_slope': dict(otype=np.ndarray, atype=np.floating,
                                   descr='Fitted slope (nslits)'),
                 'fit_b': dict(otype=np.ndarray, atype=np.floating,
                               descr='Fitted b value(nslits)'),
                 'fit_los': dict(otype=np.ndarray, atype=np.floating,
                                 descr='Fitted line width(nslits)'),
                 'resid_sky': dict(otype=np.ndarray, atype=np.floating,
                                   descr='Residuals of flexure model on sky lines (nslits)'),
                 'objra': dict(otype=np.ndarray, atype=np.floating, descr='Object RA (nslits)'),
                 'objdec': dict(otype=np.ndarray, atype=np.floating, descr='Object DEC (nslits)'),
                 'maskdef_id': dict(otype=np.ndarray, atype=np.integer, descr='Mask ID (nslits)'),
                 'rms_arc': dict(otype=np.ndarray, atype=np.floating,
                                 descr='RMS of fit (ndet, nslits)')}

    internals = ['flex_par',        # Parameters (FlexurePar)
                 'spectrograph',    # spectrograph
                 'specobjs',        # Specobjs object
                 'sobj_idx',        # (ndet, nslits); Index to specobjs (tuple of arrays)
                 'sky_table',       # Sky line table
                 # 2D models
                 'pmodel_m',
                 'pmodel_b',
                 'pmodel_l'
                ]

    def __init__(self, s1dfile=None, PYP_SPEC=None, nslits=None, det=None, 
                 SN=None, slitid=None, mn_wv=None, fit_slope=None, fit_b=None,
                 fit_los=None, objra=None, objdec=None, maskdef_id=None, rms_arc=None, 
                 resid_sky=None, indiv_fit_slope=None, indiv_fit_b=None,
                 indiv_fit_los=None):

        # Setup the DataContainer
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = {k: values[k] for k in args[1:]}
        # Init
        super().__init__(d=_d)

        # Load up specobjs
        self.specobjs = specobjs.SpecObjs.from_fitsfile(self.s1dfile, chk_version=False) 
        #  Sky lines -- This one is ASCII, so don't use load_sky_spectrum()
        sky_file = 'sky_single_mg.dat'
        self.sky_table = ascii.read(data.Paths.sky_spec / sky_file)

    # NOTE: If you make changes to how this object is bundled into the output
    # datamodel, make sure you update the documentation in
    # doc/calibrations/flexure.rst!
    def _bundle(self):
        """
        Override the base class method simply to set the HDU extension name.
        """
        return super()._bundle(ext='FLEXURE')

    def init(self, spectrograph, par):
        """ Initialize this and that about the slits, par, spectrograph
        e.g. RA, DEC, S/N

        Args:
            spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
                The spectrograph instance that sets the instrument used to take
                the observations.  Used to set :attr:`spectrograph`.
            par (:class:`~pypeit.par.pypeitpar.FlexurePar`):
                The parameters used for the flexure processing
        """
        # Internals
        self.spectrograph = spectrograph
        self.flex_par = par
        # Set
        self.PYP_SPEC = self.spectrograph.name
        self.sobj_idx = self.spectrograph.spec1d_match_spectra(self.specobjs)
        #
        self.nslits = len(self.sobj_idx[0])
        self.ndet = len(self.sobj_idx)
        
        # Fill in 1D
        self['slitid'] = self.specobjs[self.sobj_idx[0]]['SLITID'].astype(float)
        self['objra'] = self.specobjs[self.sobj_idx[0]]['RA']
        self['objdec'] = self.specobjs[self.sobj_idx[0]]['DEC']
        #self['slitname'] = self.specobjs[self.sobj_idx[0]]['MASKDEF_OBJNAME']
        self['maskdef_id'] = self.specobjs[self.sobj_idx[0]]['MASKDEF_ID']

        # Compile the list of detector *names* once
        DETs = self.specobjs.DET
        # Find which ones are actually mosaics
        is_msc = np.array([Mosaic.name_prefix in d for d in DETs]).astype(np.uint16)
        # Use the relevant parser to get the integer identifier
        det_msc_num = np.array([Mosaic.parse_name(d) if m else DetectorContainer.parse_name(d) 
                                    for d,m in zip(DETs, is_msc)])
        # Then assign the attributes
        self.is_msc = np.vstack(tuple(is_msc[self.sobj_idx[det]] for det in range(self.ndet)))
        self.det = np.vstack(tuple(det_msc_num[self.sobj_idx[det]] for det in range(self.ndet)))

        # S/N and mn_wv from the spectra
        self['SN'] = np.zeros((self.ndet, self.nslits), dtype=float)
        self['mn_wv'] = np.zeros((self.ndet, self.nslits), dtype=float)
        for det in range(self.ndet):
            self['SN'][det] = [sobj.S2N for sobj in self.specobjs[self.sobj_idx[det]]]
            self['mn_wv'][det] = [sobj.mnx_wave[0] for sobj in self.specobjs[self.sobj_idx[det]]]

    def fit_mask_surfaces(self):
        """
        Fit 2D model to linear flexure models from each slit as a function of
        RA, DEC.
        """
        # Cut on S/N
        good_SN = self['SN'] > self.flex_par['multi_min_SN']
        good_slit = np.sum(good_SN, axis=0) == self.ndet

        # Basic stats
        mu = np.median(self['indiv_fit_slope'][good_slit])
        sd = np.std(self['indiv_fit_slope'][good_slit])
        mu2 = np.median(self['indiv_fit_b'][good_slit])
        sd2 = np.std(self['indiv_fit_b'][good_slit])

        # Cut down to +/- 2sigma
        mgood = (np.abs(self['indiv_fit_slope']-mu) < 2.*sd) \
                    & ( np.abs(self['indiv_fit_b']-mu2) < 2.*sd2) & good_slit

        # Fit me (without additional rejection)
        # TODO -- Allow for x,y position instead of RA, DEC
        self.pmodel_m = fitting.robust_fit(self['objra'][mgood],
                                           self['indiv_fit_slope'][mgood], (2,2),
                                           function='polynomial2d',
                                           x2=self['objdec'][mgood])
        self.pmodel_b = fitting.robust_fit(self['objra'][mgood],
                                           self['indiv_fit_b'][mgood], (2,2),
                                           function='polynomial2d',
                                           x2=self['objdec'][mgood])
        self.pmodel_l = fitting.robust_fit(self['objra'][mgood],
                                           self['indiv_fit_los'][mgood], (2,2),
                                           function='polynomial2d',
                                           x2=self['objdec'][mgood])

    def measure_sky_lines(self):
        """Main method to analyze the sky lines for all the slits
        """

        # Init
        for key in ['indiv_fit_slope', 'indiv_fit_b', 'indiv_fit_los']:
            self[key] = np.zeros(self.nslits)

        # Loop on slits
        for i in np.arange(0,self.nslits,1):
            if (i % 10) == 0:
                msgs.info("Working on slit {} of {}".format(i, self.nslits))

            if not np.all(self['SN'][:,i] > 1.):
                continue

            # Loop on detectors
            sky_lines, sky_diffs, sky_ediffs, sky_loss = [], [], [], []
            for det in range(self.ndet):
                sobj = self.specobjs[self.sobj_idx[det][i]]

                # Measure em
                # The following will break if only boxcar...
                # TODO -- Allow for boxcar
                sky_line, sky_diff, sky_ediff, los, _ = sky_em_residuals(
                    sobj['OPT_WAVE'], 
                    sobj['OPT_COUNTS_SKY'], 
                    sobj['OPT_COUNTS_IVAR'],
                    self.sky_table['Wave'])

                # Hold em
                sky_lines.append(sky_line)
                sky_diffs.append(sky_diff)
                sky_ediffs.append(sky_ediff)
                sky_loss.append(los)

            # Concatenate
            sky_lines = np.concatenate(sky_lines)
            sky_diffs = np.concatenate(sky_diffs)
            sky_ediffs = np.concatenate(sky_ediffs)
            sky_loss = np.concatenate(sky_loss)
            
            # FIT SINGLE SLIT SKY LINES WITH A LINE           
            linear_fit = fitting.robust_fit(sky_lines,
                                            sky_diffs,
                                            weights=1./sky_ediffs**2,  
                                            function='polynomial', 
                                            order=1,
                                            maxrej=1,  # Might increase
                                            lower=3., upper=3.)
            # Save 
            self['indiv_fit_b'][i]     = linear_fit.fitc[0]
            self['indiv_fit_slope'][i] = linear_fit.fitc[1]
            self['indiv_fit_los'][i]   = np.median(sky_loss)

    def update_fit(self):
        """Update fits for each slit based on 2D model
        """
        # Do it
        self['fit_slope'] = self.pmodel_m.eval(self['objra'],x2=self['objdec'])
        self['fit_b']     = self.pmodel_b.eval(self['objra'],x2=self['objdec'])
        self['fit_los']   = self.pmodel_l.eval(self['objra'],x2=self['objdec'])

        # CALCULATE RESIDUALS FROM FIT
        #   Only for QA (I think)
        resid_sky = []
        for i in range(self.nslits):

            # Require sufficient S/N in reddest detector
            if self['SN'][-1,i] > 0:
                # Load up the full spectrum
                tmp_wave, all_flux, all_sky, all_ivar = np.ndarray(0), \
                    np.ndarray(0), np.ndarray(0), np.ndarray(0)
                # TODO -- Allow for Boxcar
                for det in range(self.ndet):
                    sobj = self.specobjs[self.sobj_idx[det][i]]
                    tmp_wave = np.concatenate((tmp_wave, sobj.OPT_WAVE))
                    all_flux = np.concatenate((all_flux, sobj.OPT_COUNTS))
                    all_sky = np.concatenate((all_sky, sobj.OPT_COUNTS_SKY))
                    all_ivar = np.concatenate((all_ivar, sobj.OPT_COUNTS_IVAR))
                
                # Massage
                fitwave  = self['fit_slope'][i]*tmp_wave + self['fit_b'][i]
                all_wave = tmp_wave - fitwave

                # TRIM ENDS
                all_wave=all_wave[5:-15]
                all_flux=all_flux[5:-15]
                all_ivar=all_ivar[5:-15]
                all_sky=all_sky[5:-15]

                # REMOVE CRAZY 500-SIGMA VALUES
                cmask = (all_sky > np.percentile(all_sky,0.1)) & (all_sky < np.percentile(all_sky,99.9))

                m=np.median(all_sky[cmask])
                s=np.std(all_sky[cmask])
                mm = (all_sky > 500.*s + m) | (all_sky < m-50.*s)
                all_sky[mm] = m
                all_ivar[mm] = 1e6
                if (np.sum(mm) > 10):
                    msgs.warn('Removing more than 10 pixels of data')
                
                _,diff,diff_err,_,_ = sky_em_residuals(all_wave, all_sky, all_ivar,
                                                       self.sky_table['Wave'])
                m = np.isfinite(diff)
                sky_mean = np.average(np.abs(diff[m]), weights = 1./diff_err[m]**2)
                resid_sky = np.append(resid_sky,sky_mean)

            else:
                resid_sky = np.append(resid_sky,-1)

        self['resid_sky'] = resid_sky

    def qa_plots(self, plot_dir:str, root:str):
        """Generate QA plots

        Args:
            plot_dir (str): Top-lvel folder for QA
                QA/ is generated beneath this, as needed
            root (str): Root for output files
        """

        # Generate QA folder as need be
        qa_dir = pathlib.Path(plot_dir) / 'QA'
        if not qa_dir.is_dir():
            qa_dir.mkdir(parents=True)
        
        '''
        # Slopes
        pdf2 = matplotlib.backends.backend_pdf.PdfPages(os.path.join(qa_dir, 'flex_slits_'+root+'.pdf'))
        plt.rcParams.update({'figure.max_open_warning': 0})
        for i in np.arange(0,self.nslits,1):

            if not np.all(self['SN'][:,i] > 0.):
                continue


            # SKY LINES FIRST
            r_sky_line, r_sky_diff,r_sky_ediff,r_los,r_elos = sky_em_residuals(hdu[r].data['OPT_WAVE'], \
                                                    hdu[r].data['OPT_COUNTS_SKY'],\
                                                    hdu[r].data['OPT_COUNTS_IVAR'])

            b_sky_line, b_sky_diff,b_sky_ediff,b_los,b_elos = sky_em_residuals(hdu[b].data['OPT_WAVE'], \
                                                    hdu[b].data['OPT_COUNTS_SKY'],\
                                                    hdu[b].data['OPT_COUNTS_IVAR'])

            fig, (ax1,ax2) = plt.subplots(1, 2,figsize=(20,4))
            ax1.plot(r_sky_line,r_sky_diff,'ro',alpha=0.8,label='Red chip: Sky Emission')
            ax1.plot(b_sky_line,b_sky_diff,'bo',alpha=0.8,label='Blue chip: Sky Emission')
            ax1.errorbar(b_sky_line,b_sky_diff,yerr=b_sky_ediff,fmt='none',ecolor='b',alpha=0.5)
            ax1.errorbar(r_sky_line,r_sky_diff,yerr=r_sky_ediff,fmt='none',ecolor='r',alpha=0.5)
            ax1.text(6320,0,'{}'.format(b),fontsize=11)
            ax1.text(8500,0,'{}'.format(r),fontsize=11)
            ax1.set_ylim(-0.45,0.45)

            x=np.arange(6000,9000,1)
            l1 = slits['fit_slope'][i]*x + slits['fit_b'][i]
            l2 = fslits['fit_slope'][i]*x + fslits['fit_b'][i]
            ax1.plot(x,l1,'-')
            ax1.plot(x,l2,'--')
            ax1.axhline(linewidth=1, color='grey',alpha=0.5)
            ax1.set_ylabel('Wavelength offset (AA)')
            ax1.set_xlabel('Wavelength (AA)')
            ax1.set_xlim(6300,9100)
            t = 'Sky Line Fits , resid = {:0.4f} AA, arc = {:0.2f}'.format(slits['resid_sky'][i],0.32*slits['rms_arc_r'][i])
            ax1.set_title(t)

            sky_diff  = np.concatenate((r_sky_diff,b_sky_diff),axis=None)
            sky_lines = np.concatenate((r_sky_line,b_sky_line),axis=None)
            sky_ediff = np.concatenate((r_sky_ediff,b_sky_ediff),axis=None)
            sky_los   = np.concatenate((r_los,b_los),axis=None)


            ax2.plot(r_sky_line,r_los,'ro',alpha=0.8,label='Red chip: Sky Emission')
            ax2.plot(b_sky_line,b_los,'bo',alpha=0.8,label='Blue chip: Sky Emission')
            ax2.errorbar(r_sky_line,r_los,yerr=r_elos,fmt='none',ecolor='r',alpha=0.5)
            ax2.errorbar(b_sky_line,b_los,yerr=b_elos,fmt='none',ecolor='b',alpha=0.5)
            ax2.axhline(fslits['fit_los'][i],linewidth=1, color='grey',alpha=0.5)

            ax2.set_title('Line widths')
            ax2.set_xlabel('Wavelength (AA)')
            ax2.set_ylim(0.3,0.8)
            ax2.set_xlim(6300,9100)

            pdf2.savefig()
        pdf2.close()
        plt.close('all')
        '''

        #########################################################################
        # CREATE FULL MASK FITS
        pdf = matplotlib.backends.backend_pdf.PdfPages(
            plot_dir+'QA/flex_mask_'+root+'.pdf')
        xslit = self['objra']
        yslit = self['objdec']
        t=2.

        mu =  np.median(self['indiv_fit_slope'])
        sd =  np.std(self['indiv_fit_slope'])
        mu2 =  np.median(self['indiv_fit_b'])
        sd2 =  np.std(self['indiv_fit_b'])
        mu3 =  np.median(self['indiv_fit_los'])
        sd3 =  np.std(self['indiv_fit_los'])

        # PLOT FITTED VALUES
        fig, (ax1,ax2,ax3) = plt.subplots(1, 3,figsize=(22,5))
    
        mm1=-0.00005
        mm2=0.00005
        print(mu-t*sd,mu+t*sd)
        ax1.scatter(xslit,yslit,c=self['indiv_fit_slope'],
                    cmap="cool",vmin = mm1,vmax=mm2 )# mu-t*sd,vmax=mu+t*sd)
        ax1.set_ylabel('Dec [deg]')
        ax1.set_xlabel('RA [deg]')
        ax1.set_title('Wave MEASURE: line slope')
        #cax, _ = matplotlib.colorbar.make_axes(ax1)
        #normalize = matplotlib.colors.Normalize(vmin = mu-t*sd,vmax=mu+t*sd)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='cool',norm=normalize)


        ax2.scatter(xslit,yslit,c=self['indiv_fit_b'],cmap="summer",
                    vmin = mu2-t*sd2,vmax=mu2+t*sd2)
        ax2.set_ylabel('Dec [deg]')
        ax2.set_xlabel('RA [deg]')
        ax2.set_title('Wave MEASURE: line intercept')
        cax, _ = matplotlib.colorbar.make_axes(ax2)
        normalize = matplotlib.colors.Normalize(vmin = mu2-t*sd2,vmax=mu2+t*sd2)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='summer',norm=normalize)


        ax3.scatter(xslit,yslit,c=self['indiv_fit_los'],cmap="cool",vmin = mu3-t*sd3,vmax=mu3+t*sd3)
        ax3.set_ylabel('Dec [deg]')
        ax3.set_xlabel('RA [deg]')
        ax3.set_title('Wave MEASURE: line width')
        cax, _ = matplotlib.colorbar.make_axes(ax3)
        normalize = matplotlib.colors.Normalize(vmin = mu3-t*sd3,vmax=mu3+t*sd3)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='cool',norm=normalize)

        pdf.savefig()
        
        #######################
        # PLOT MEASURED VALUES
        fig, (ax1,ax2,ax3) = plt.subplots(1, 3,figsize=(22,5))
    
        ax1.scatter(xslit,yslit,c=self['fit_slope'],
                    cmap="cool",vmin = mu-t*sd,vmax=mu+t*sd)

        ax1.set_ylabel('Dec [deg]')
        ax1.set_xlabel('RA [deg]')
        ax1.set_title('Wave fit: line slope')
        cax, _ = matplotlib.colorbar.make_axes(ax1)
        normalize = matplotlib.colors.Normalize(vmin = mu-t*sd,vmax=mu+t*sd)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='cool',norm=normalize)


        ax2.scatter(xslit,yslit,c=self['fit_b'],
                    cmap="summer",vmin = mu2-t*sd2,vmax=mu2+t*sd2)
        ax2.set_ylabel('Dec [deg]')
        ax2.set_xlabel('RA [deg]')
        ax2.set_title('Wave fit: line intercept')
        cax, _ = matplotlib.colorbar.make_axes(ax2)
        normalize = matplotlib.colors.Normalize(vmin = mu2-t*sd2,vmax=mu2+t*sd2)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='summer',norm=normalize)


        ax3.scatter(xslit,yslit,c=self['fit_los'],
                    cmap="cool",vmin = mu3-t*sd3,vmax=mu3+t*sd3)
        ax3.set_ylabel('Dec [deg]')
        ax3.set_xlabel('RA [deg]')
        ax3.set_title('Wave fit: line width')
        cax, _ = matplotlib.colorbar.make_axes(ax3)
        normalize = matplotlib.colors.Normalize(vmin = mu3-t*sd3,vmax=mu3+t*sd3)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='cool',norm=normalize)

        
        pdf.close()
