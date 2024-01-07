#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the viewing of a processed FITS file
with extras.  Run above the Science/ folder.
"""

import os
import argparse
import numpy as np

from astropy.table import Table

from pypeit import msgs
from pypeit import io
from pypeit import slittrace
from pypeit.core import gui
from pypeit.core.parse import get_dnum


def parse_args(options=None, return_parser=False):

    parser = argparse.ArgumentParser(description='Display sky subtracted, spec2d image in the'
                                                 'interactive object finding GUI.  Run above'
                                                 'the Science/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default=None, help='PYPEIT spec2d file')
    parser.add_argument("--list", default=False, help="List the extensions only?",
                        action="store_true")
    parser.add_argument('--det', default=1, type=int, help="Detector")
    parser.add_argument("--old", default=False, action="store_true", help="Used old slit tracing")

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def parse_traces(hdulist_1d, det_nm):
    """Extract the relevant trace information
    """
    traces = dict(traces=[], fwhm=[])
    pkflux = []
    for hdu in hdulist_1d:
        if det_nm in hdu.name:
            tbl = Table(hdu.data)
            trace = tbl['TRACE_SPAT']
            fwhm = tbl['FWHMFIT']
            obj_id = hdu.name.split('-')[0]
            traces['traces'].append(trace.copy())
            traces['fwhm'].append(np.median(fwhm))
            pkflux.append(np.median(tbl['BOX_COUNTS']))
    traces['pkflux'] = np.array(pkflux)
    return traces


def main(args):

    raise NotImplementedError('This script is currently out of date.')

    # List only?
    hdu = io.fits_open(args.file)
    head0 = hdu[0].header
    if args.list:
        hdu.info()
        return

    # Init
    sdet = get_dnum(args.det, prefix=False)

    # One detector, sky sub for now
    names = [hdu[i].name for i in range(len(hdu))]

    try:
        exten = names.index('DET{:s}-PROCESSED'.format(sdet))
    except:  # Backwards compatability
        msgs.error('Requested detector {:s} was not processed.\n'
                   'Maybe you chose the wrong one to view?\n'
                   'Set with --det= or check file contents with --list'.format(sdet))
    sciimg = hdu[exten].data
    try:
        exten = names.index('DET{:s}-SKY'.format(sdet))
    except:  # Backwards compatability
        msgs.error('Requested detector {:s} has no sky model.\n'
                   'Maybe you chose the wrong one to view?\n'
                   'Set with --det= or check file contents with --list'.format(sdet))
    skymodel = hdu[exten].data
    try:
        exten = names.index('DET{:s}-MASK'.format(sdet))
    except ValueError:  # Backwards compatability
        msgs.error('Requested detector {:s} has no bit mask.\n'
                   'Maybe you chose the wrong one to view?\n'
                   'Set with --det= or check file contents with --list'.format(sdet))

    mask = hdu[exten].data
    frame = (sciimg - skymodel) * (mask == 0)

    mdir = head0['PYPMFDIR']
    mkey = head0['FRAMMKEY']
    mast_key = '{0}_{1:02d}'.format(mkey, args.det)
    if not os.path.exists(mdir):
        mdir_base = os.path.join(os.getcwd(), os.path.basename(mdir))
        msgs.warn('Master file dir: {0} does not exist. Using {1}'.format(mdir, mdir_base))
        mdir = mdir_base

    # Assumes a MasterSlit file has been written
    #slits = slittrace.SlitTraceSet.from_master('{0}_{1:02d}'.format(head0['TRACMKEY'], args.det),
    #                                           mdir)
    # Load the slits information
    slits = slittrace.SlitTraceSet.from_master(mast_key, mdir)

    # Object traces
    left, right, mask = slits.select_edges()
    msgs.error("You need to choose which slits you care about here")

    # Get object traces
    spec1d_file = args.file.replace('spec2d', 'spec1d')
    if os.path.isfile(spec1d_file):
        hdulist_1d = io.fits_open(spec1d_file)
    else:
        hdulist_1d = []
        msgs.warn('Could not find spec1d file: {:s}'.format(spec1d_file) + msgs.newline() +
                  '                          No objects were extracted.')

    msgs.error("This code needs to be refactored since tslits_dict was removed...")
    import pdb
    pdb.set_trace()
    tslits_dict['objtrc'] = parse_traces(hdulist_1d, det_nm)
    obj_trace = parse_traces(hdulist_1d, 'DET{:s}'.format(sdet))

    # TODO :: Need to include standard star trace in the spec2d files
    std_trace = None

    # Extract some trace models
    fwhm = 2  # Start with some default value
    # TODO: Dictionaries like this are a pet peeve of mine.  I'd prefer
    # either individual objects or a class with a well-formed data model.
    # TODO: Why do all of these dictionary elements need fwhm?  Can they
    # be different?
    trace_models = dict()
    # Brightest object on slit
    trace_models['object'] = dict(trace_model=None, fwhm=fwhm)
    if len(obj_trace['pkflux']) > 0:
        smash_peakflux = obj_trace['pkflux']
        ibri = smash_peakflux.argmax()
        trace_models['object']['trace_model'] = obj_trace['traces'][ibri]
        trace_models['object']['fwhm'] = obj_trace['fwhm'][ibri]
    # Standard star trace
    trace_models['std'] = dict(trace_model=std_trace, fwhm=trace_models['object']['fwhm'])
    # Trace of the slit edge
    # TODO: Any particular reason to use the lefts?
    trace_models['slit'] = dict(trace_model=left.copy(), fwhm=trace_models['object']['fwhm'])

    # Finally, initialise the GUI
    gui.object_find.initialise(args.det, frame, left, right, obj_trace, trace_models, None,
                               printout=True, slit_ids=slits.id)

    ofgui = gui_object_find.initialise(args.det, frame, tslits_dict, None, printout=True, slit_ids=slits.id)



def illum_profile_spatial(self, skymask=None, trim_edg=(0, 0), debug=False):
    """
    Calculate the residual spatial illumination profile using the sky regions.

    The residual is calculated using the differential:

    .. code-block:: python

        correction = amplitude * (1 + spatial_shift * (dy/dx)/y)

    where ``y`` is the spatial profile determined from illumflat, and
    spatial_shift is the residual spatial flexure shift in units of pixels.

     Args:
        skymask (`numpy.ndarray`_):
            Mask of sky regions where the spatial illumination will be determined
        trim_edg (:obj:`tuple`):
            A tuple of two ints indicated how much of the slit edges should be
            trimmed when fitting to the spatial profile.
        debug (:obj:`bool`):
            Show debugging plots?
    """

    msgs.info("Performing spatial sensitivity correction")
    # Setup some helpful parameters
    skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)
    hist_trim = 0  # Trim the edges of the histogram to take into account edge effects
    gpm = self.sciImg.select_flag(invert=True)
    slitid_img_init = self.slits.slit_img(pad=0, initial=True, flexure=self.spat_flexure_shift)
    spatScaleImg = np.ones_like(self.sciImg.image)
    # For each slit, grab the spatial coordinates and a spline
    # representation of the spatial profile from the illumflat
    rawimg = self.sciImg.image.copy()
    numbins = int(np.max(self.slits.get_slitlengths(initial=True, median=True)))
    spatbins = np.linspace(0.0, 1.0, numbins + 1)
    spat_slit = 0.5 * (spatbins[1:] + spatbins[:-1])
    slitlength = np.median(self.slits.get_slitlengths(median=True))
    coeff_fit = np.zeros((self.slits.nslits, 2))
    for sl, slitnum in enumerate(self.slits.spat_id):
        msgs.info("Deriving spatial correction for slit {0:d}/{1:d}".format(sl + 1, self.slits.spat_id.size))
        # Get the initial slit locations
        onslit_b_init = (slitid_img_init == slitnum)

        # Synthesize ximg, and edgmask from slit boundaries. Doing this outside this
        # routine would save time. But this is pretty fast, so we just do it here to make the interface simpler.
        spatcoord, edgmask = pixels.ximg_and_edgemask(self.slits_left[:, sl], self.slits_right[:, sl],
                                                      onslit_b_init, trim_edg=trim_edg)

        # Make the model histogram
        xspl = np.linspace(0.0, 1.0, 10 * int(slitlength))  # Sub sample each pixel with 10 subpixels
        # TODO: caliBrate is no longer a dependency. If you need these b-splines pass them in.
        modspl = self.caliBrate.flatimages.illumflat_spat_bsplines[sl].value(xspl)[0]
        gradspl = interpolate.interp1d(xspl, np.gradient(modspl) / modspl, kind='linear', bounds_error=False,
                                       fill_value='extrapolate')

        # Ignore skymask
        coord_msk = onslit_b_init & gpm
        hist, _ = np.histogram(spatcoord[coord_msk], bins=spatbins, weights=rawimg[coord_msk])
        cntr, _ = np.histogram(spatcoord[coord_msk], bins=spatbins)
        hist_slit_all = hist / (cntr + (cntr == 0))
        histmod, _ = np.histogram(spatcoord[coord_msk], bins=spatbins, weights=gradspl(spatcoord[coord_msk]))
        hist_model = histmod / (cntr + (cntr == 0))

        # Repeat with skymask
        coord_msk = onslit_b_init & gpm & skymask_now
        hist, _ = np.histogram(spatcoord[coord_msk], bins=spatbins, weights=rawimg[coord_msk])
        cntr, _ = np.histogram(spatcoord[coord_msk], bins=spatbins)
        hist_slit = hist / (cntr + (cntr == 0))

        # Prepare for fit - take the non-zero elements and trim slit edges
        if hist_trim == 0:
            ww = (hist_slit != 0)
            xfit = spat_slit[ww]
            yfit = hist_slit_all[ww]
            mfit = hist_model[ww]
        else:
            ww = (hist_slit[hist_trim:-hist_trim] != 0)
            xfit = spat_slit[hist_trim:-hist_trim][ww]
            yfit = hist_slit_all[hist_trim:-hist_trim][ww]
            mfit = hist_model[hist_trim:-hist_trim][ww]

        # Fit the function
        spat_func = lambda par, ydata, model: par[0]*(1 + par[1] * model) - ydata
        res_lsq = least_squares(spat_func, [np.median(yfit), 0.0], args=(yfit, mfit))
        spatnorm = spat_func(res_lsq.x, 0.0, gradspl(spatcoord[onslit_b_init]))
        spatnorm /= spat_func(res_lsq.x, 0.0, gradspl(0.5))
        # Set the scaling factor
        spatScaleImg[onslit_b_init] = spatnorm
        coeff_fit[sl, :] = res_lsq.x

    if debug:
        from matplotlib import pyplot as plt
        xplt = np.arange(24)
        plt.subplot(121)
        plt.plot(xplt[0::2], coeff_fit[::2, 0], 'rx')
        plt.plot(xplt[1::2], coeff_fit[1::2, 0], 'bx')
        plt.subplot(122)
        plt.plot(xplt[0::2], coeff_fit[::2, 1]/10, 'rx')
        plt.plot(xplt[1::2], coeff_fit[1::2, 1]/10, 'bx')
        plt.show()
        plt.imshow(spatScaleImg, vmin=0.99, vmax=1.01)
        plt.show()
        plt.subplot(133)
        plt.plot(xplt[0::2], coeff_fit[::2, 2], 'rx')
        plt.plot(xplt[1::2], coeff_fit[1::2, 2], 'bx')
        plt.show()
    # Apply the relative scale correction
    self.apply_relative_scale(spatScaleImg)

def illum_profile_spectral(self, global_sky, skymask=None):
    """Calculate the residual spectral illumination profile using the sky regions.
    This uses the same routine as the flatfield spectral illumination profile.

     Args:
         global_sky (`numpy.ndarray`_):
            Model of the sky
         skymask (`numpy.ndarray`_, optional):
            Mask of sky regions where the spectral illumination will be determined
    """
    trim = self.par['calibrations']['flatfield']['slit_trim']
    sl_ref = self.par['calibrations']['flatfield']['slit_illum_ref_idx']
    smooth_npix = self.par['calibrations']['flatfield']['slit_illum_smooth_npix']
    gpm = self.sciImg.select_flag(invert=True)
    # Note :: Need to provide wavelength to illum_profile_spectral (not the tilts) so that the
    # relative spectral sensitivity is calculated at a given wavelength for all slits simultaneously.
    scaleImg = flatfield.illum_profile_spectral(self.sciImg.image.copy(), self.waveimg, self.slits,
                                                slit_illum_ref_idx=sl_ref, model=global_sky, gpmask=gpm,
                                                skymask=skymask, trim=trim, flexure=self.spat_flexure_shift,
                                                smooth_npix=smooth_npix)
    # Now apply the correction to the science frame
    self.apply_relative_scale(scaleImg)



def global_skysub_old(self, skymask=None, update_crmask=True, trim_edg=(0,0),
                  previous_sky=None, show_fit=False, show=False, show_objs=False, objs_not_masked=False,
                  reinit_bpm:bool=True):
    """
    Perform global sky subtraction. This SlicerIFU-specific routine ensures that the
    edges of the slits are not trimmed, and performs a spatial and spectral
    correction using the sky spectrum, if requested. See Reduce.global_skysub()
    for parameter definitions.

    See base class method for description of parameters.

    Args:
        reinit_bpm (:obj:`bool`, optional):
            If True (default), the bpm is reinitialized to the initial bpm
            Should be False on the final run in case there was a failure
            upstream and no sources were found in the slit/order
    """
    if self.par['reduce']['findobj']['skip_skysub']:
        msgs.info("Skipping global sky sub as per user request")
        return np.zeros_like(self.sciImg.image)

    # Generate a global sky sub for all slits separately
    global_sky_sep = super().global_skysub(skymask=skymask, update_crmask=update_crmask,
                                           trim_edg=trim_edg, show_fit=show_fit, show=show,
                                           show_objs=show_objs, reinit_bpm=reinit_bpm)
    # Check if any slits failed
    if np.any(global_sky_sep[self.slitmask >= 0] == 0) and not self.bkg_redux:
        # Cannot continue without a sky model for all slits
        msgs.error("Global sky subtraction has failed for at least one slit.")

    # Check if flexure or a joint fit is requested
    if not self.par['reduce']['skysub']['joint_fit'] and self.par['flexure']['spec_method'] == 'skip':
        return global_sky_sep
    if self.wv_calib is None:
        msgs.error("A wavelength calibration is needed (wv_calib) if a joint sky fit is requested.")
    msgs.info("Generating wavelength image")

    self.waveimg = self.wv_calib.build_waveimg(self.tilts, self.slits, spat_flexure=self.spat_flexure_shift)
    # Calculate spectral flexure
    method = self.par['flexure']['spec_method']
    # TODO :: Perhaps include a new label for IFU flexure correction - e.g. 'slitcen_relative' or 'slitcenIFU' or 'IFU'
    #      :: If a new label is introduced, change the other instances of 'method' (see below), and in flexure.spec_flexure_qa()
    if method in ['slitcen']:
        self.slitshift = self.calculate_flexure(global_sky_sep)
        # Recalculate the wavelength image, and the global sky taking into account the spectral flexure
        msgs.info("Generating wavelength image, now accounting for spectral flexure")
        self.waveimg = self.wv_calib.build_waveimg(self.tilts, self.slits, spec_flexure=self.slitshift,
                                                   spat_flexure=self.spat_flexure_shift)


    # If the joint fit or spec/spat sensitivity corrections are not being performed, return the separate slits sky
    if not self.par['reduce']['skysub']['joint_fit']:
        return global_sky_sep

    # Do the spatial scaling first
    # if self.par['scienceframe']['process']['use_illumflat']:
    #     # Perform the correction
    #     self.illum_profile_spatial(skymask=skymask)
    #     # Re-generate a global sky sub for all slits separately
    #     global_sky_sep = Reduce.global_skysub(self, skymask=skymask, update_crmask=update_crmask, trim_edg=trim_edg,
    #                                           show_fit=show_fit, show=show, show_objs=show_objs)

    # Use sky information in all slits to perform a joint sky fit
    global_sky = self.joint_skysub(skymask=skymask, update_crmask=update_crmask, trim_edg=trim_edg,
                                   show_fit=show_fit, show=show, show_objs=show_objs,
                                   objs_not_masked=objs_not_masked)

    return global_sky
