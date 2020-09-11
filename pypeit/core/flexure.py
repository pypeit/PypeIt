""" Module for flexure routines

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import inspect

import numpy as np
import copy
from matplotlib import pyplot as plt
from matplotlib import gridspec

from astropy import stats
from astropy import units
import scipy.signal
import scipy.optimize as opt
from scipy import interpolate

from linetools.spectra import xspectrum1d

from pypeit import msgs
from pypeit import utils
from pypeit.display import display
from pypeit.core import arc
from pypeit.core import qa
from pypeit.core import fitting
from pypeit.core.moment import moment1d

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
        if debug:
            embed(header='68 of flexure')
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
        embed(header='83 of flexure.py')

    return lag_max[0]



def load_sky_spectrum(sky_file):
    """
    Load a sky spectrum into an XSpectrum1D object

    ..todo -- Try to eliminate the XSpectrum1D dependancy

    Args:
        sky_file: str

    Returns:
        sky_spec: XSpectrum1D
          spectrum
    """
    return xspectrum1d.XSpectrum1D.from_file(sky_file)


def spec_flex_shift(obj_skyspec, arx_skyspec, arx_lines, mxshft=20):
    """ Calculate shift between object sky spectrum and archive sky spectrum

    Args:
        obj_skyspec (:class:`linetools.spectra.xspectrum1d.XSpectrum1d`):
            Spectrum of the sky related to our object
        arx_skyspec (:class:`linetools.spectra.xspectrum1d.XSpectrum1d`):
            Archived sky spectrum
        arx_lines (tuple): Line information returned by arc.detect_lines for
            the Archived sky spectrum
        mxshft (float, optional):
            Maximum allowed shift from flexure;  note there are cases that
            have been known to exceed even 30 pixels..

    Returns:
        dict: Contains flexure info
    """

    # TODO None of these routines should have dependencies on XSpectrum1d!

    # Determine the brightest emission lines
    msgs.warn("If we use Paranal, cut down on wavelength early on")
    arx_amp, arx_amp_cont, arx_cent, arx_wid, _, arx_w, arx_yprep, nsig \
            = arx_lines
    obj_amp, obj_amp_cont, obj_cent, obj_wid, _, obj_w, obj_yprep, nsig_obj \
            = arc.detect_lines(obj_skyspec.flux.value)

    # Keep only 5 brightest amplitude lines (xxx_keep is array of
    # indices within arx_w of the 5 brightest)
    arx_keep = np.argsort(arx_amp[arx_w])[-5:]
    obj_keep = np.argsort(obj_amp[obj_w])[-5:]

    # Calculate wavelength (Angstrom per pixel)
    arx_disp = np.append(arx_skyspec.wavelength.value[1]-arx_skyspec.wavelength.value[0],
                         arx_skyspec.wavelength.value[1:]-arx_skyspec.wavelength.value[:-1])
    obj_disp = np.append(obj_skyspec.wavelength.value[1]-obj_skyspec.wavelength.value[0],
                         obj_skyspec.wavelength.value[1:]-obj_skyspec.wavelength.value[:-1])

    # Calculate resolution (lambda/delta lambda_FWHM)..maybe don't need
    # this? can just use sigmas
    arx_idx = (arx_cent+0.5).astype(np.int)[arx_w][arx_keep]   # The +0.5 is for rounding
    arx_res = arx_skyspec.wavelength.value[arx_idx]/\
              (arx_disp[arx_idx]*(2*np.sqrt(2*np.log(2)))*arx_wid[arx_w][arx_keep])
    obj_idx = (obj_cent+0.5).astype(np.int)[obj_w][obj_keep]   # The +0.5 is for rounding
    obj_res = obj_skyspec.wavelength.value[obj_idx]/ \
              (obj_disp[obj_idx]*(2*np.sqrt(2*np.log(2)))*obj_wid[obj_w][obj_keep])

    if not np.all(np.isfinite(obj_res)):
        msgs.warn('Failed to measure the resolution of the object spectrum, likely due to error '
                   'in the wavelength image.')
        return None
    msgs.info("Resolution of Archive={0} and Observation={1}".format(np.median(arx_res),
                                                                     np.median(obj_res)))

    # Determine sigma of gaussian for smoothing
    arx_sig2 = np.power(arx_disp[arx_idx]*arx_wid[arx_w][arx_keep], 2)
    obj_sig2 = np.power(obj_disp[obj_idx]*obj_wid[obj_w][obj_keep], 2)

    arx_med_sig2 = np.median(arx_sig2)
    obj_med_sig2 = np.median(obj_sig2)

    if obj_med_sig2 >= arx_med_sig2:
        smooth_sig = np.sqrt(obj_med_sig2-arx_med_sig2)  # Ang
        smooth_sig_pix = smooth_sig / np.median(arx_disp[arx_idx])
        arx_skyspec = arx_skyspec.gauss_smooth(smooth_sig_pix*2*np.sqrt(2*np.log(2)))
    else:
        msgs.warn("Prefer archival sky spectrum to have higher resolution")
        smooth_sig_pix = 0.
        msgs.warn("New Sky has higher resolution than Archive.  Not smoothing")
        #smooth_sig = np.sqrt(arx_med_sig**2-obj_med_sig**2)

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
    if np.any(np.isfinite(corr[subpix_grid.astype(np.int)])):
        fit = fitting.PypeItFit(xval=subpix_grid, yval=corr[subpix_grid.astype(np.int)],
                                func='polynomial', order=np.atleast_1d(2))
        fit.fit()
        success = True
        max_fit = -0.5 * fit.fitc[1] / fit.fitc[2]
    else:
        fit = fitting.PypeItFit(xval=subpix_grid, yval=0.0*subpix_grid,
                                funct='polynomial', order=np.atleast_1d(2))
        fit.fit()
        success = False
        max_fit = 0.0
        msgs.warn('Flexure compensation failed for one of your objects')

    #Calculate and apply shift in wavelength
    shift = float(max_fit)-lag0
    msgs.info("Flexure correction of {:g} pixels".format(shift))
    #model = (fit[2]*(subpix_grid**2.))+(fit[1]*subpix_grid)+fit[0]

    return dict(polyfit=fit, shift=shift, subpix=subpix_grid,
                corr=corr[subpix_grid.astype(np.int)], sky_spec=obj_skyspec, arx_spec=arx_skyspec,
                corr_cen=corr.size/2, smooth=smooth_sig_pix, success=success)


def flexure_interp(sky_wave, fdict, xin=None):
    """
    Apply interpolation with the flexure dict

    Args:
        sky_wave (`numpy.ndarray`_):
            Wavelengths of the extracted sky
        fdict (dict):
            Holds the various flexure items
        xin (`numpy.ndarray`_):
            spectral values (in pixels) to apply flexure correction (normalized to a value between 0-1)

    Returns:
        xspectrum1d.XSpectrum1D:  New sky spectrum (mainly for QA)

    """
    # Simple interpolation to apply
    npix = len(sky_wave)
    x = np.linspace(0., 1., npix)
    # Apply
    msgs.info("Applying flexure correction")
    f = interpolate.interp1d(x, sky_wave, bounds_error=False, fill_value="extrapolate")
    if xin is None: return f(x + fdict['shift'] / (npix - 1))
    else: return f(xin + fdict['shift'] / (npix - 1))


def spec_flexure_slit(waveimg, skyimg, traces, slits, slitmask, bpm, sky_file, mxshft=None):
    """Correct wavelengths for flexure, slit by slit

    Args:
        waveimg (`numpy.ndarray`_):
            A 2D array providing the wavelength (in A) of each pixel of an image
        skyimg (`numpy.ndarray`_):
            A 2D array of the global sky model for each slit
        traces (`numpy.ndarray`_):
            The central trace of each slit
        slits (:class:`pypeit.slittrace.SlitTraceSet`_):
            Slit trace set
        slitmask (`numpy.ndarray`_):
            Slit mask
        bpm (`numpy.ndarray`_):
            True = masked slit
        sky_file (str):
            Sky file
        mxshft (int, optional):
            Passed to flex_shift()

    Returns:
        numpy.ndarray: The input waveimg, corrected for spectral flexure
        list:  list of dicts containing flexure results of each slit
            Filled with a basically empty dict if the slit is skipped
    """
    msgs.work("Consider doing 2 passes in flexure as in LowRedux")
    waveimg_flex_corr = waveimg.copy()

    # Load Archive. Save the line information to avoid the performance hit from calling it on the archive sky spectrum
    # multiple times
    sky_spectrum = load_sky_spectrum(sky_file)
    sky_lines = arc.detect_lines(sky_spectrum.flux.value)

    nslits = len(bpm)
    gpm = np.invert(bpm)
    gdslits = np.where(gpm)[0]

    # Loop on objects
    flex_list = []

    # Loop over slits, and then over objects here
    nspec = slits.nspec
    spec_vec = np.arange(nspec)
    for islit in range(nslits):
        msgs.info("Working on spectral flexure of slit: {:d}".format(islit))

        # Reset
        flex_dict = dict(polyfit=[], shift=[], subpix=[], corr=[],
                         corr_cen=[], spec_file=sky_file, smooth=[],
                         arx_spec=[], sky_spec=[])
        # If no objects on this slit append an empty dictionary
        if islit not in gdslits:
            flex_list.append(flex_dict.copy())
            continue

        slit_spat = slits.spat_id[islit]
        mask = (slitmask == slit_spat) & gpm

        # Fill in the boxcar extraction tags
        box_denom = moment1d(waveimg * mask > 0.0, traces[:, islit], 2, row=spec_vec)[0]
        wghts = (box_denom + (box_denom == 0.0))
        sky_flux = moment1d(skyimg * mask, traces[:, islit], 2, row=spec_vec)[0] / wghts
        # Denom is computed in case the trace goes off the edge of the image
        sky_wave = moment1d(waveimg * mask, traces[:, islit], 2, row=spec_vec)[0] / wghts

        # Generate 1D spectrum for object
        obj_sky = xspectrum1d.XSpectrum1D.from_tuple((sky_wave, sky_flux))

        # Calculate the shift
        fdict = spec_flex_shift(obj_sky, sky_spectrum, sky_lines, mxshft=mxshft)
        # Failed?
        if fdict is not None:
            # Update dict
            for key in ['polyfit', 'shift', 'subpix', 'corr', 'corr_cen', 'smooth', 'arx_spec']:
                flex_dict[key].append(fdict[key])
            # Interpolate
            flex_dict['sky_spec'].append(xspectrum1d.XSpectrum1D.from_tuple((flexure_interp(sky_wave, fdict), sky_flux)))
            # Now correct the waveimg
            waveimg_nrm = (waveimg[slitmask == slit_spat] - sky_wave.min()) / (sky_wave.max() - sky_wave.min())
            waveimg_flex_corr[slitmask == slit_spat] = flexure_interp(sky_wave, fdict, xin=waveimg_nrm)

        # Append, this will be an empty dictionary if the flexure failed
        flex_list.append(flex_dict.copy())

    return waveimg_flex_corr, flex_list


def spec_flexure_corrQA(ax, this_flex_dict, cntr, name):
    # Fit
    fit = this_flex_dict['polyfit'][cntr]
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


def spec_flexure_qa(slitords, bpm, basename, det, flex_list, specobjs=None, out_dir=None):
    """

    Args:
        slitords (`numpy.ndarray`_):
            Array of slit/order numbers
        bpm (`numpy.ndarray`_):
            True = masked slit
        basename (str):
        det (int):
        flex_list (list):
        specobjs:
        out_dir:

    """
    plt.rcdefaults()
    plt.rcParams['font.family'] = 'times new roman'

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

        # Parse and Setup
        if slit_cen:
            nobj = 1
            ncol = 1
        else:
            indx = specobjs.slitorder_indices(slitord)
            this_specobjs = specobjs[indx]
            nobj = np.sum(indx)
            if nobj == 0: continue
            ncol = min(3, nobj)
        this_flex_dict = flex_list[islit]

        # Check that the default was overwritten
        if len(this_flex_dict['shift']) == 0:
            continue

        nrow = nobj // ncol + ((nobj % ncol) > 0)
        # Outfile, one QA file per slit
        outfile = qa.set_qa_filename(basename, method + '_corr', det=det, slit=slitord, out_dir=out_dir)
        plt.figure(figsize=(8, 5.0))
        plt.clf()
        gs = gridspec.GridSpec(nrow, ncol)
        # TODO -- This cntr is crummy and needs to be replaced by a DataContainer
        #  for flex_dict and flex_list
        cntr = 0
        # Correlation QA
        if slit_cen:
            ax = plt.subplot(gs[0, 0])
            spec_flexure_corrQA(ax, this_flex_dict, cntr, 'Slit Center')
        else:
            for specobj in this_specobjs:
                if specobj is None or (specobj.BOX_WAVE is None and specobj.OPT_WAVE is None):
                    continue
                ax = plt.subplot(gs[cntr//ncol, cntr % ncol])
                spec_flexure_corrQA(ax, this_flex_dict, cntr, '{:s}'.format(specobj.NAME))
                cntr += 1
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
        outfile = qa.set_qa_filename(basename, method+'_sky', det=det,slit=slitord, out_dir=out_dir)
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
        msgs.info("Wrote spectral flexure QA: {}".format(outfile))
        #plt.close()

    plt.rcdefaults()


def calculate_image_offset(image, im_ref, nfit=3):
    """Calculate the x,y offset between two images

    Args:
        image (`numpy.ndarray`_):
            Image that we want to measure the shift of (relative to im_ref)
        im_ref (`numpy.ndarray`_):
            Reference image
        nfit (int, optional):
            Number of pixels (left and right of the maximum) to include in
            fitting the peak of the cross correlation.

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
    # Subtract median (should be close to zero, anyway)
    image -= np.median(image)
    im_ref -= np.median(im_ref)

    # cross correlate (note, convolving seems faster)
    ccorr = scipy.signal.correlate2d(im_ref, image, boundary='fill', mode='same')
    #ccorr = scipy.signal.fftconvolve(im_ref, image[::-1, ::-1], mode='same')

    # Find the maximum
    amax = np.unravel_index(np.argmax(ccorr), ccorr.shape)

    # Perform a 2D Gaussian fit
    x = np.arange(amax[0]-nfit, amax[0] + nfit+1)
    y = np.arange(amax[1]-nfit, amax[1] + nfit+1)
    initial_guess = (np.max(ccorr), amax[0], amax[1], 3, 3, 0, 0)
    xx, yy = np.meshgrid(x, y, indexing='ij')

    # Fit the neighborhood of the maximum to calculate the offset
    popt, _ = opt.curve_fit(fitting.twoD_Gaussian, (xx, yy),
                            ccorr[amax[0]-nfit:amax[0]+nfit+1, amax[1]-nfit:amax[1]+nfit+1].ravel(),
                            p0=initial_guess)
    # Return the RA and DEC shift, in pixels
    return popt[1] - ccorr.shape[0]//2, popt[2] - ccorr.shape[1]//2

