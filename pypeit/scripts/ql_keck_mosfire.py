#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script runs PypeIt on a pair of NIRES images (A-B)
"""

import os
import argparse
import numpy as np
import copy

from astropy.io import fits
from astropy.table import Table
from pypeit import pypeit
from pypeit import par, msgs
from pypeit import pypeitsetup
from pypeit import wavecalib
from pypeit import wavetilts
from pypeit import spec2dobj
from pypeit import coadd2d
from pypeit import specobjs
from pypeit import slittrace
from pypeit import reduce
from pypeit import calibrations
from pypeit.display import display
from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.parse import get_dnum
from pypeit.core.wavecal import wvutils
from pypeit import sensfunc
from pypeit.core import flux_calib
from astropy.stats import sigma_clipped_stats
from IPython import embed


# A trick from stackoverflow to allow multi-line output in the help:
#https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

def parse_args(options=None, return_parser=False):

    parser = argparse.ArgumentParser(description='Script to run PypeIt on a pair of MOSFIRE files (A-B)', formatter_class=SmartFormatter)
    parser.add_argument('full_rawpath', type=str, help='Full path to the raw files')
    parser.add_argument('files', type=str, nargs='+', help='list of frames i.e. img1.fits img2.fits')
    parser.add_argument('--samp_fact', default=1.0, type=float,
                        help="Make the wavelength grid finer (samp_fact > 1.0) or coarser (samp_fact < 1.0) by this sampling factor")
    parser.add_argument("--flux", default=False, action='store_true',
                        help="This option will multiply in sensitivity function to obtain a flux calibrated 2d spectrum")
    parser.add_argument("--mask_cr", default=False, action='store_true',
                        help="This option turns on cosmic ray rejection. This improves the reduction but doubles runtime.")
    parser.add_argument('--box_radius', type=float, help='Set the radius for the boxcar extraction')
    parser.add_argument('--offset', type=float, default=None,
                        help='R|Override the automatic offsets determined from the headers. Offset is in pixels.\n'
                        'This option is useful if a standard dither pattern was not executed.\n'
                        'The offset convention is such that a negative offset will move the (negative) B image to the left')
    parser.add_argument("--redux_path", type=str, default=os.getcwd(),
                        help="Location where reduction outputs should be stored.")
    parser.add_argument("--master_dir", type=str, default=os.getenv('MOSFIRE_MASTERS'),
                        help="Location of PypeIt Master files used for the reduction.")
    parser.add_argument('--embed', default=False, help='Upon completion embed in ipython shell',
                        action='store_true')
    parser.add_argument("--show", default=False, action="store_true",
                        help="Show the reduction steps. Equivalent to the -s option when running pypeit.")

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def config_lines(args):

    # Config the run
    cfg_lines = ['[rdx]']
    cfg_lines += ['    spectrograph = {0}'.format('keck_mosfire')]
    cfg_lines += ['    redux_path = {0}'.format(args.redux_path)]
    cfg_lines += ['    scidir = Science_QL']
    # Calibrations
    cfg_lines += ['[baseprocess]']
    cfg_lines += ['    use_pixelflat = False']
    cfg_lines += ['    use_illumflat = False']
    cfg_lines += ['[calibrations]']
    cfg_lines += ['    [[wavelengths]]']
    cfg_lines += ['        refframe = observed']
    if not args.mask_cr:
        cfg_lines += ['[scienceframe]']
        cfg_lines += ['    [[process]]']
        cfg_lines += ['        mask_cr = False']
    cfg_lines += ['[reduce]']
    cfg_lines += ['    [[extraction]]']
    cfg_lines += ['        skip_optimal = True']
    if args.box_radius is not None: # Boxcar radius
        cfg_lines += ['        boxcar_radius = {0}'.format(args.box_radius)]
    cfg_lines += ['    [[findobj]]']
    cfg_lines += ['        skip_second_find = True']

    return cfg_lines


def save_exposure(fitstbl, frame, spectrograph, science_path, par, caliBrate, all_spec2d, all_specobjs):
    """
    Save the outputs from extraction for a given exposure

    Args:
        frame (:obj:`int`):
            0-indexed row in the metadata table with the frame
            that has been reduced.
        all_spec2d(:class:`pypeit.spec2dobj.AllSpec2DObj`):
        sci_dict (:obj:`dict`):
            Dictionary containing the primary outputs of
            extraction
        basename (:obj:`str`):
            The root name for the output file.

    Returns:
        None or SpecObjs:  All of the objects saved to disk

    """
    # TODO: Need some checks here that the exposure has been reduced?

    # Get the basename
    basename = fitstbl.construct_basename(frame)

    # Determine the headers
    row_fitstbl = fitstbl[frame]
    # Need raw file header information
    rawfile = fitstbl.frame_paths(frame)
    head2d = fits.getheader(rawfile, ext=spectrograph.primary_hdrext)

    # Check for the directory
    if not os.path.isdir(science_path):
        os.makedirs(science_path)

    subheader = spectrograph.subheader_for_spec(row_fitstbl, head2d)
    # 1D spectra
    if all_specobjs.nobj > 0:
        # Spectra
        outfile1d = os.path.join(science_path, 'spec1d_{:s}.fits'.format(basename))
        all_specobjs.write_to_fits(subheader, outfile1d,
                                   update_det=par['rdx']['detnum'],
                                   slitspatnum=par['rdx']['slitspatnum'])
        # Info
        outfiletxt = os.path.join(science_path, 'spec1d_{:s}.txt'.format(basename))
        all_specobjs.write_info(outfiletxt, spectrograph.pypeline)
    else:
        outfile1d = None

    # 2D spectra
    outfile2d = os.path.join(science_path, 'spec2d_{:s}.fits'.format(basename))
    # Build header
    pri_hdr = all_spec2d.build_primary_hdr(head2d, spectrograph,
                                           redux_path=par['rdx']['redux_path'],
                                           master_key_dict=caliBrate.master_key_dict,
                                           master_dir=caliBrate.master_dir,
                                           subheader=subheader)
    # Write
    all_spec2d.write_to_fits(outfile2d, pri_hdr=pri_hdr, update_det=par['rdx']['detnum'])
    return outfile2d, outfile1d

def run_pair(A_files, B_files, caliBrate, spectrograph, det, parset, show=False, std_trace=None):
    """
    Peform 2d extraction for a set of files at the same unique A-B offset location.

    Parameters
    ----------
    A_files (list of strings):
       Files at A position for this offset
    B_files (list of strings)
       Files at B position for this offeset
    caliBrate (object):
       CaliBrate object
    spectrograph (object):
       spectrograph object
    det (int):
       Detector number
    parset (parsect object)
       Parset
    show (bool, optional):
       Show 2d reduction outputs. Default=False
    std_trace (string, optional)
       Trace for standard star. Default=None

    Returns
    -------
    spec2DObj_A, spec2DObj_B

    spec2DObj_A (object, Spec2D):
       Spec2d Object for extraction at A position
    spec2DObj_B (object, Spec2D)
       Spec2d Object for extraction at B position

    """

    # Build Science image
    sciImg = buildimage.buildimage_fromlist(
        spectrograph, det, parset['scienceframe'], list(A_files), bpm=caliBrate.msbpm, slits=caliBrate.slits, ignore_saturation=False)

    # Background Image?
    sciImg = sciImg.sub(buildimage.buildimage_fromlist(
        spectrograph, det, parset['scienceframe'], list(B_files), bpm=caliBrate.msbpm, slits=caliBrate.slits, ignore_saturation=False),
        parset['scienceframe']['process'])
    # Instantiate Reduce object
    # Required for pypeline specific object
    # At instantiaton, the fullmask in self.sciImg is modified
    redux = reduce.Reduce.get_instance(sciImg, spectrograph, parset, caliBrate, 'science', ir_redux=True, show=show,
                                       det=det)

    skymodel, objmodel, ivarmodel, outmask, sobjs, scaleimg, waveimg, tilts = redux.run(
        std_trace=std_trace, return_negative=True, show_peaks=show)

    # TODO -- Do this upstream
    # Tack on detector
    for sobj in sobjs:
        sobj.DETECTOR = sciImg.detector

    # Construct table of spectral flexure
    spec_flex_table = Table()
    spec_flex_table['spat_id'] = caliBrate.slits.spat_id
    spec_flex_table['sci_spec_flexure'] = redux.slitshift

    # Construct the Spec2DObj with the positive image
    spec2DObj_A = spec2dobj.Spec2DObj(det=det,
                                      sciimg=sciImg.image,
                                      ivarraw=sciImg.ivar,
                                      skymodel=skymodel,
                                      objmodel=objmodel,
                                      ivarmodel=ivarmodel,
                                      scaleimg=scaleimg,
                                      waveimg=waveimg,
                                      bpmmask=outmask,
                                      detector=sciImg.detector,
                                      sci_spat_flexure=sciImg.spat_flexure,
                                      sci_spec_flexure=spec_flex_table,
                                      vel_corr=None,
                                      vel_type=parset['calibrations']['wavelengths']['refframe'],
                                      tilts=tilts,
                                      slits=copy.deepcopy(caliBrate.slits))
    spec2DObj_A.process_steps = sciImg.process_steps
    all_spec2d = spec2dobj.AllSpec2DObj()
    all_spec2d['meta']['ir_redux'] = True
    all_spec2d[det] = spec2DObj_A
    # Save image A but with all the objects extracted, i.e. positive and negative
    # outfile2d, outfile1d = save_exposure(fitstbl, 0, spectrograph, science_path, parset, caliBrate, all_spec2d, sobjs)

    # Construct the Spec2DObj with the negative image
    spec2DObj_B = spec2dobj.Spec2DObj(det=det,
                                      sciimg=-sciImg.image,
                                      ivarraw=sciImg.ivar,
                                      skymodel=-skymodel,
                                      objmodel=-objmodel,
                                      ivarmodel=ivarmodel,
                                      scaleimg=scaleimg,
                                      waveimg=waveimg,
                                      bpmmask=outmask,
                                      detector=sciImg.detector,
                                      sci_spat_flexure=sciImg.spat_flexure,
                                      sci_spec_flexure=spec_flex_table,
                                      vel_corr=None,
                                      vel_type=parset['calibrations']['wavelengths']['refframe'],
                                      tilts=tilts,
                                      slits=copy.deepcopy(caliBrate.slits))
    return spec2DObj_A, spec2DObj_B

def main(args):


    # Calibration Master directory
    if args.master_dir is None:
        msgs.error('You need to set an environment variable MOSFIRE_MASTERS that points at the Master Calibs')

    # Define some hard wired master files here to be later parsed out of the directory
    slit_masterframe_name = os.path.join(args.master_dir, 'MasterSlits_D_8191_01.fits.gz')
    tilts_masterframe_name = os.path.join(args.master_dir, 'MasterTilts_D_1_01.fits')
    wvcalib_masterframe_name = os.path.join(args.master_dir, 'MasterWaveCalib_D_1_01.fits')
    std_spec1d_file = os.path.join(args.master_dir, 'spec1d_m201024_0232-gd71_MOSFIRE_2020Oct24T153823.555.fits')
    sensfunc_masterframe_name = os.path.join(args.master_dir, 'sens_m201024_0232-gd71_MOSFIRE_2020Oct24T153823.555.fits')
    if (not os.path.isfile(slit_masterframe_name) or  not os.path.isfile(tilts_masterframe_name) or \
        not os.path.isfile(tilts_masterframe_name) or not os.path.isfile(sensfunc_masterframe_name) or \
        not os.path.isfile(std_spec1d_file)):
        msgs.error('Master frames not found. Check that environment variable MOSFIRE_MASTERS  points at the Master Calibs')

    # Read in the spectrograph, config the parset
    spectrograph = load_spectrograph('keck_mosfire')
    spectrograph_def_par = spectrograph.default_pypeit_par()
    parset = par.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                          merge_with=config_lines(args))
    science_path = os.path.join(parset['rdx']['redux_path'], parset['rdx']['scidir'])

    # Parse the files sort by MJD
    files = np.array([os.path.join(args.full_rawpath, file) for file in args.files])
    nfiles = len(files)
    mjds = np.zeros(nfiles)
    for ifile, file in enumerate(files):
        hdr = fits.getheader(file, spectrograph.primary_hdrext)
        try:
           mjds[ifile] = hdr['MJD-OBS']
        except:
            msgs.warn('File {:} has no MJD in header'.format(file))
            mjds[ifile] = 0
    files = files[np.argsort(mjds)]

    # We need the platescale
    platescale = spectrograph.get_detector_par(None, 1)['platescale']
    # Parse the offset information out of the headers. TODO in the future get this out of fitstable
    dither_pattern, dither_id, offset_arcsec = spectrograph.parse_dither_pattern(files)
    if len(np.unique(dither_pattern)) > 1:
        msgs.error('Currently this script is supported only for a single type of dither pattern')
    A_files = files[dither_id == 'A']
    B_files = files[dither_id == 'B']
    nA = len(A_files)
    nB = len(B_files)

    # Print out a report on the offsets
    msg_string = msgs.newline()  +     '****************************************************'
    msg_string += msgs.newline() +     ' Summary of offsets for dither pattern:   {:s}'.format(dither_pattern[0])
    msg_string += msgs.newline() +     '****************************************************'
    msg_string += msgs.newline() +     'filename     Position         arcsec    pixels    '
    msg_string += msgs.newline() +     '----------------------------------------------------'
    for iexp, file in enumerate(files):
        msg_string += msgs.newline() + '    {:s}    {:s}   {:6.2f}    {:6.2f}'.format(
            os.path.basename(file), dither_id[iexp], offset_arcsec[iexp], offset_arcsec[iexp]/platescale)
    msg_string += msgs.newline() +     '****************************************************'
    msgs.info(msg_string)

    #offset_dith_pix = offset_dith_pix = offset_arcsec_A[0]/sciImg.detector.platescale

    ## Read in the master frames that we need
    ##
    det = 1 # MOSFIRE has a single detector
    if std_spec1d_file is not None:
        # Get the standard trace if need be
        sobjs = specobjs.SpecObjs.from_fitsfile(std_spec1d_file)
        this_det = sobjs.DET == det
        if np.any(this_det):
            sobjs_det = sobjs[this_det]
            sobjs_std = sobjs_det.get_std()
            std_trace = None if sobjs_std is None else sobjs_std.TRACE_SPAT.flatten()
        else:
            std_trace = None
    else:
        std_trace = None

    # Read in the msbpm
    sdet = get_dnum(det, prefix=False)
    msbpm = spectrograph.bpm(A_files[0], det)
    # Read in the slits
    slits = slittrace.SlitTraceSet.from_file(slit_masterframe_name)
    # Reset the bitmask
    slits.mask = slits.mask_init.copy()
    # Read in the wv_calib
    wv_calib = wavecalib.WaveCalib.from_file(wvcalib_masterframe_name)
    #wv_calib.is_synced(slits)
    slits.mask_wvcalib(wv_calib)
    # Read in the tilts
    tilts_obj = wavetilts.WaveTilts.from_file(tilts_masterframe_name)
    tilts_obj.is_synced(slits)
    slits.mask_wavetilts(tilts_obj)

    # Build the Calibrate object
    caliBrate = calibrations.Calibrations(None, parset['calibrations'], spectrograph, None)
    caliBrate.slits = slits
    caliBrate.msbpm = msbpm
    caliBrate.wavetilts = tilts_obj
    caliBrate.wv_calib = wv_calib

    # Find the unique throw absolute value, which defines each MASK_NOD seqeunce
    #uniq_offsets, _ = np.unique(offset_arcsec, return_inverse=True)
    uniq_throws, uni_indx = np.unique(np.abs(offset_arcsec), return_inverse=True)
    # uniq_throws = uniq values of the dither throw
    # uni_indx = indices into the uniq_throws array needed to reconstruct the original array
    nuniq = uniq_throws.size
    spec2d_list =[]
    offset_ref = offset_arcsec[0]
    offsets_dith_pix = []
    # Generalize to a multiple slits, doing one slit at a time?
    islit = 0
    # Loop over the unique throws and create a spec2d_A and spec2D_B for each, which are then
    # fed into coadd2d with the correct offsets
    for iuniq in range(nuniq):
        A_ind = (uni_indx == iuniq) & (dither_id == 'A')
        B_ind = (uni_indx == iuniq) & (dither_id == 'B')
        A_files = files[A_ind]
        B_files = files[B_ind]
        A_offset = offset_arcsec[A_ind]
        B_offset = offset_arcsec[B_ind]
        throw = np.abs(A_offset[0])
        msgs.info('Reducing A-B pairs for throw = {:}'.format(throw))
        spec2DObj_A, spec2DObj_B = run_pair(A_files, B_files, caliBrate, spectrograph, det, parset,
                                            show=args.show, std_trace=std_trace)
        spec2d_list += [spec2DObj_A, spec2DObj_B]
        offsets_dith_pix += [(np.mean(A_offset) - offset_ref)/platescale, (np.mean(B_offset) - offset_ref)/platescale]

    offsets_dith_pix = np.array(offsets_dith_pix)
    #else:
    #    msgs.error('Unrecognized mode')

    if args.offset is not None:
        offsets_pixels = np.array([0.0, args.offset])
        msgs.info('Using user specified offsets instead: {:5.2f}'.format(args.offset))
    else:
        offsets_pixels = offsets_dith_pix

    # Instantiate Coadd2d
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_list, spectrograph, parset, det=det,
                                         offsets=offsets_pixels, weights='uniform', ir_redux=True,
                                         debug=args.show, samp_fact=args.samp_fact)
    # Coadd the slits
    coadd_dict_list = coadd.coadd(only_slits=None, interp_dspat=False)  # TODO implement only_slits later
    # Create the pseudo images
    pseudo_dict = coadd.create_pseudo_image(coadd_dict_list)

    if args.flux:
        # Load the sensitivity function
        wave_sens, sfunc, _, _, _ = sensfunc.SensFunc.load(sensfunc_masterframe_name)
        # Interpolate the sensitivity function onto the wavelength grid of the data. Since the image is rectified
        # this is trivial and we don't need to do a 2d interpolation
        sens_factor = flux_calib.get_sensfunc_factor(
            pseudo_dict['wave_mid'][:,islit], wave_sens, sfunc, fits.getheader(files[0])['TRUITIME'],
            extrap_sens=parset['fluxcalib']['extrap_sens'])
        # Compute the median sensitivity and set the sensitivity to zero at locations 100 times the median. This
        # prevents the 2d image from blowing up where the sens_factor explodes because there is no throughput
        sens_gpm = sens_factor < 100.0*np.median(sens_factor)
        sens_factor_masked = sens_factor*sens_gpm
        sens_factor_img = np.repeat(sens_factor_masked[:, np.newaxis], pseudo_dict['nspat'], axis=1)
        imgminsky = sens_factor_img*pseudo_dict['imgminsky']
        imgminsky_gpm = sens_gpm[:, np.newaxis] & pseudo_dict['inmask']
    else:
        imgminsky= pseudo_dict['imgminsky']

    ##########################
    # Now display the images #
    ##########################
    display.connect_to_ginga(raise_err=True, allow_new=True)
    # TODO: Bug in ginga prevents me from using cuts here for some reason
    mean, med, sigma = sigma_clipped_stats(imgminsky[imgminsky_gpm], sigma_lower=3.0, sigma_upper=3.0)
    chname_skysub = 'fluxed-skysub-det{:s}'.format(sdet) if args.flux else 'skysub-det{:s}'.format(sdet)
    cuts_skysub = (med - 3.0 * sigma, med + 3.0 * sigma)
    cuts_resid = (-5.0, 5.0)
    #fits.writeto('/Users/joe/ginga_test.fits',imgminsky, overwrite=True)
    #fits.writeto('/Users/joe/ginga_mask.fits',imgminsky_gpm.astype(float), overwrite=True)
    #embed()


    # Clear all channels at the beginning
    # TODO: JFH For some reason Ginga crashes when I try to put cuts in here.
    viewer, ch_skysub = display.show_image(imgminsky, chname=chname_skysub, waveimg=pseudo_dict['waveimg'],
                                   clear=True, cuts= cuts_skysub)
    slit_left, slit_righ, _ = pseudo_dict['slits'].select_edges()
    slit_id = slits.slitord_id[0]
    display.show_slits(viewer, ch_skysub, slit_left, slit_righ, slit_ids=slit_id)

    # SKRESIDS
    chname_skyresids = 'sky_resid-det{:s}'.format(sdet)
    image = pseudo_dict['imgminsky']*np.sqrt(pseudo_dict['sciivar']) * pseudo_dict['inmask']  # sky residual map
    viewer, ch_skyresids = display.show_image(image, chname_skyresids, waveimg=pseudo_dict['waveimg'],
                                  cuts=cuts_resid)

    display.show_slits(viewer, ch_skyresids, slit_left, slit_righ, slit_ids=slits.slitord_id[0])
    shell = viewer.shell()
    out = shell.start_global_plugin('WCSMatch')
    out = shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [chname_skysub], {})


    # TODO extract along a spatial position



    if args.embed:
        embed()

    return 0

