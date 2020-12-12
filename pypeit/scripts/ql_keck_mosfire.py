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
from pypeit import display
from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.parse import get_dnum

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
    parser.add_argument('-A','--Afiles', type=str, nargs='+', help='list of frames at dither position A, i.e. -A A1.fits A2.fits')
    parser.add_argument('-B','--Bfiles', type=str, nargs='+', help='list of frames at dither position B  i.e. -B B1.fits B2.fits')
    parser.add_argument('--samp_fact', default=1.0, type=float,
                        help="Make the wavelength grid finer (samp_fact > 1.0) or coarser (samp_fact < 1.0) by this sampling factor")
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
    cfg_lines += ['        frame = observed']
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

def parse_dither_pattern(file_list, ext):

    nfiles = len(file_list)
    offset_arcsec = np.zeros(nfiles)
    dither_pattern = []
    dither_id = []
    for ifile, file in enumerate(file_list):
        hdr = fits.getheader(file, ext)
        dither_pattern.append(hdr['PATTERN'])
        dither_id.append(hdr['FRAMEID'])
        offset_arcsec[ifile] = hdr['YOFFSET']

    return dither_pattern, dither_id, offset_arcsec



def main(args):


    # Build the fitstable since we currently need it for output. This should not be the case!
    A_files = [os.path.join(args.full_rawpath, file) for file in args.Afiles]
    B_files = [os.path.join(args.full_rawpath, file) for file in args.Bfiles]
    data_files = A_files + B_files
    ps = pypeitsetup.PypeItSetup(A_files, path='./', spectrograph_name='keck_mosfire')
    ps.build_fitstbl()
    fitstbl = ps.fitstbl

    # Read in the spectrograph, config the parset
    spectrograph = load_spectrograph('keck_mosfire')
    spectrograph_def_par = spectrograph.default_pypeit_par()
    parset = par.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                          merge_with=config_lines(args))
    science_path = os.path.join(parset['rdx']['redux_path'], parset['rdx']['scidir'])

    # Calibration Master directory
    if args.master_dir is None:
        msgs.error("You need to set an Environmental variable MOSFIRE_MASTERS that points at the Master Calibs")

    # Define some hard wired master files here to be later parsed out of the directory
    slit_masterframe_name = os.path.join(args.master_dir, 'MasterSlits_E_15_01.fits.gz')
    tilts_masterframe_name = os.path.join(args.master_dir, 'MasterTilts_E_1_01.fits')
    wvcalib_masterframe_name = os.path.join(args.master_dir, 'MasterWaveCalib_E_1_01.fits')
    # For now don't require a standard
    std_outfile=None
    #std_outfile = os.path.join('/Users/joe/Dropbox/PypeIt_Redux/MOSFIRE/Nov19/quicklook/Science/',
    #                           'spec1d_m191118_0064-GD71_MOSFIRE_2019Nov18T104704.507.fits')
    # make the get_std from pypeit a utility function or class method
    det = 1 # MOSFIRE has a single detector
    if std_outfile is not None:
        # Get the standard trace if need be
        sobjs = specobjs.SpecObjs.from_fitsfile(std_outfile)
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
    wv_calib.is_synced(slits)
    slits.mask_wvcalib(wv_calib)
    # Read in the tilts
    tilts_obj = wavetilts.WaveTilts.from_file(tilts_masterframe_name)
    tilts_obj.is_synced(slits)
    slits.mask_wavetilts(tilts_obj)


    # Build Science image
    sciImg = buildimage.buildimage_fromlist(spectrograph, det, parset['scienceframe'], A_files, bpm=msbpm,
                                            slits=slits, ignore_saturation=False)

    # Background Image?
    sciImg = sciImg.sub(buildimage.buildimage_fromlist(spectrograph, det,  parset['scienceframe'], B_files, bpm=msbpm,
                                                       slits=slits,ignore_saturation=False), parset['scienceframe']['process'])
    # Build the Calibrate object
    caliBrate = calibrations.Calibrations(None, parset['calibrations'], spectrograph, None)
    caliBrate.slits = slits
    caliBrate.wavetilts = tilts_obj
    caliBrate.wv_calib = wv_calib

    # Instantiate Reduce object
    # Required for pypeline specific object
    # At instantiaton, the fullmask in self.sciImg is modified
    redux = reduce.Reduce.get_instance(sciImg, spectrograph, parset, caliBrate, 'science', ir_redux=True, show=args.show,
                                       det=det, std_outfile=std_outfile)

    manual_extract_dict = None
    skymodel, objmodel, ivarmodel, outmask, sobjs, waveImg, tilts = redux.run(
        std_trace=std_trace, return_negative=True, manual_extract_dict=manual_extract_dict, show_peaks=args.show)

    # TODO -- Do this upstream
    # Tack on detector
    for sobj in sobjs:
        sobj.DETECTOR = sciImg.detector

    # Construct the Spec2DObj with the positive image
    spec2DObj_A = spec2dobj.Spec2DObj(det=det,
                                      sciimg=sciImg.image,
                                      ivarraw=sciImg.ivar,
                                      skymodel=skymodel,
                                      objmodel=objmodel,
                                      ivarmodel=ivarmodel,
                                      waveimg=waveImg,
                                      bpmmask=outmask,
                                      detector=sciImg.detector,
                                      sci_spat_flexure=sciImg.spat_flexure,
                                      tilts=tilts,
                                      slits=copy.deepcopy(caliBrate.slits))
    spec2DObj_A.process_steps = sciImg.process_steps
    all_spec2d = spec2dobj.AllSpec2DObj()
    all_spec2d['meta']['ir_redux'] = True
    all_spec2d[det] = spec2DObj_A
    # Save image A but with all the objects extracted, i.e. positive and negative
    #outfile2d, outfile1d = save_exposure(fitstbl, 0, spectrograph, science_path, parset, caliBrate, all_spec2d, sobjs)

    # Construct the Spec2DObj with the negative image
    spec2DObj_B= spec2dobj.Spec2DObj(det=det,
                                     sciimg=-sciImg.image,
                                     ivarraw=sciImg.ivar,
                                     skymodel=-skymodel,
                                     objmodel=-objmodel,
                                     ivarmodel=ivarmodel,
                                     waveimg=waveImg,
                                     bpmmask=outmask,
                                     detector=sciImg.detector,
                                     sci_spat_flexure=sciImg.spat_flexure,
                                     tilts=tilts,
                                     slits=copy.deepcopy(caliBrate.slits))

    # Parse the offset information out of the headers. TODO in the future get this out of fitstable
    dither_pattern_A, dither_id_A, offset_arcsec_A = parse_dither_pattern(A_files, spectrograph.primary_hdrext)
    dither_pattern_B, dither_id_B, offset_arcsec_B = parse_dither_pattern(B_files, spectrograph.primary_hdrext)
    # Print out a report on the offsets
    msg_string = msgs.newline()  +     '****************************************************'
    msg_string += msgs.newline() +     ' Summary of offsets for dither pattern:   {:s}'.format(dither_pattern_A[0])
    msg_string += msgs.newline() +     '****************************************************'
    msg_string += msgs.newline() +     'Position     filename         arcsec    pixels    '
    msg_string += msgs.newline() +     '----------------------------------------------------'
    for iexp, file in enumerate(A_files):
        msg_string += msgs.newline() + '    A    {:s}   {:6.2f}    {:6.2f}'.format(
            os.path.basename(file), offset_arcsec_A[iexp], offset_arcsec_A[iexp]/sciImg.detector.platescale)
    for iexp, file in enumerate(B_files):
        msg_string += msgs.newline() + '    B    {:s}   {:6.2f}    {:6.2f}'.format(
            os.path.basename(file), offset_arcsec_B[iexp], offset_arcsec_B[iexp]/sciImg.detector.platescale)
    msg_string += msgs.newline() +     '****************************************************'
    msgs.info(msg_string)

    #offset_dith_pix = offset_dith_pix = offset_arcsec_A[0]/sciImg.detector.platescale
    offsets_dith_pix = (np.array([0.0,np.mean(offset_arcsec_B) - np.mean(offset_arcsec_A)]))/sciImg.detector.platescale
    if args.offset is not None:
        offsets_pixels = np.array([0.0, args.offset])
        msgs.info('Using user specified offsets instead: {:5.2f}'.format(args.offset))
    else:
        offsets_pixels = offsets_dith_pix


    spec2d_list = [spec2DObj_A, spec2DObj_B]
    # Instantiate Coadd2d
    coadd = coadd2d.CoAdd2D.get_instance(spec2d_list, spectrograph, parset, det=det,
                                         offsets=offsets_pixels, weights='uniform', ir_redux=True,
                                         debug=args.show, samp_fact=args.samp_fact)
    # Coadd the slits
    coadd_dict_list = coadd.coadd(only_slits=None, interp_dspat=False)  # TODO implement only_slits later
    # Create the pseudo images
    pseudo_dict = coadd.create_pseudo_image(coadd_dict_list)

    ##########################
    # Now display the images #
    ##########################
    display.display.connect_to_ginga(raise_err=True, allow_new=True)
    # Bug in ginga prevents me from using cuts here for some reason
    #mean, med, sigma = sigma_clipped_stats(pseudo_dict['imgminsky'][pseudo_dict['inmask']], sigma_lower=5.0,sigma_upper=5.0)
    #cut_min = mean - 4.0 * sigma
    #cut_max = mean + 4.0 * sigma
    chname_skysub='skysub-det{:s}'.format(sdet)
    # Clear all channels at the beginning
    # TODO: JFH For some reason Ginga crashes when I try to put cuts in here.
    viewer, ch = ginga.show_image(pseudo_dict['imgminsky'], chname=chname_skysub, waveimg=pseudo_dict['waveimg'],
                                   clear=True) # cuts=(cut_min, cut_max),
    slit_left, slit_righ, _ = pseudo_dict['slits'].select_edges()
    slit_id = slits.slitord_id[0]
    ginga.show_slits(viewer, ch, slit_left, slit_righ, slit_ids=slit_id)

    # SKRESIDS
    chname_skyresids = 'sky_resid-det{:s}'.format(sdet)
    image = pseudo_dict['imgminsky']*np.sqrt(pseudo_dict['sciivar']) * pseudo_dict['inmask']  # sky residual map
    viewer, ch = ginga.show_image(image, chname_skyresids, waveimg=pseudo_dict['waveimg'],
                                  cuts=(-5.0, 5.0),)
    ginga.show_slits(viewer, ch, slit_left, slit_righ, slit_ids=slits.slitord_id[0])
    shell = viewer.shell()
    out = shell.start_global_plugin('WCSMatch')
    out = shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [chname_skyresids], {})

    if args.embed:
        embed()

    return 0

