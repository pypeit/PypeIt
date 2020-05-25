#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script runs PypeIt on a pair of NIRES images (A-B)
"""
import argparse

from pypeit import msgs
import os
import sys
import numpy as np
import copy
from astropy.io import fits
from IPython import embed

from pypeit import pypeit
from pypeit import par, msgs
from pypeit import pypeitsetup
from pypeit import wavecalib
from pypeit import wavetilts
from pypeit import spec2dobj
from pypeit.core import framematch
from pypeit import specobjs
from pypeit.images import buildimage
from pypeit import slittrace
from pypeit.spectrographs.util import load_spectrograph
from pypeit import reduce
from pypeit import calibrations
import warnings

def parser(options=None):

    parser = argparse.ArgumentParser(description='Script to run PypeIt on a pair of MOSFIRE files (A-B)')
    parser.add_argument('full_rawpath', type=str, help='Full path to the raw files')
    parser.add_argument('fileA', type=str, help='A frame')
    parser.add_argument('fileB', type=str, help='B frame')
    parser.add_argument('-b', '--box_radius', type=float, help='Set the radius for the boxcar extraction')
    parser.add_argument("--redux_path", type=str, default=os.getcwd(),
                        help="Location where reduction outputs should be stored.")
    parser.add_argument("--show", default=False, action="store_true",
                        help="Show the reduction steps. Equivalent to the -s option when running pypeit.")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs

def config_lines(pargs):

    # Config the run
    cfg_lines = ['[rdx]']
    cfg_lines += ['    spectrograph = {0}'.format('keck_mosfire')]
    cfg_lines += ['    redux_path = {0}'.format(pargs.redux_path)]
    cfg_lines += ['    scidir = Science_QL']
    # Calibrations
    cfg_lines += ['[baseprocess]']
    cfg_lines += ['    use_pixelflat = False']
    cfg_lines += ['    use_illumflat = False']
    cfg_lines += ['[calibrations]']
    cfg_lines += ['    [[wavelengths]]']
    cfg_lines += ['        frame = observed']
    cfg_lines += ['[scienceframe]']
    cfg_lines += ['    [[process]]']
    cfg_lines += ['        mask_cr = False']
    cfg_lines += ['[reduce]']
    cfg_lines += ['    [[extraction]]']
    cfg_lines += ['        skip_optimal = True']
    if pargs.box_radius is not None: # Boxcar radius
        cfg_lines += ['        boxcar_radius = {0}'.format(pargs.box_radius)]
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



def main(pargs):

    # Build the fitstable since we currently need it for output. This should not be the case!
    data_files = [os.path.join(pargs.full_rawpath, pargs.fileA), os.path.join(pargs.full_rawpath,pargs.fileB)]
    ps = pypeitsetup.PypeItSetup([data_files[0]], path='./', spectrograph_name='keck_mosfire')
    ps.build_fitstbl()
    fitstbl = ps.fitstbl

    # Read in the spectrograph, config the parset
    spectrograph = load_spectrograph('keck_mosfire')
    spectrograph_def_par = spectrograph.default_pypeit_par()
    parset = par.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                          merge_with=config_lines(pargs))
    science_path = os.path.join(parset['rdx']['redux_path'], parset['rdx']['scidir'])

    # Calibration Master directory
    master_dir = os.getenv('MOSFIRE_MASTERS')
    if master_dir is None:
        msgs.error("You need to set an Environmental variable MOSFIRE_MASTERS that points at the Master Calibs")

    # Define some hard wired master files here to be later parsed out of the directory
    slit_masterframe_name = os.path.join(master_dir, 'MasterSlits_A_3_01.fits.gz')
    tilts_masterframe_name = os.path.join(master_dir, 'MasterTilts_A_2_01.fits')
    wvcalib_masterframe_name = os.path.join(master_dir, 'MasterWaveCalib_A_2_01.fits')
    std_outfile =os.path.join('/Users/joe/Dropbox/PypeIt_Redux/MOSFIRE/Nov19/quicklook/Science/',
                              'spec1d_m191118_0064-GD71_MOSFIRE_2019Nov18T104704.507.fits')

    sci_files = [os.path.join(pargs.full_rawpath, pargs.fileA)]
    bg_files  = [os.path.join(pargs.full_rawpath, pargs.fileB)]

    # Read in the msbpm
    det = 1 # MOSFIRE has a single detector
    msbpm = spectrograph.bpm(sci_files[0], det)
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

    # Get the standard trace if need be
    sobjs = specobjs.SpecObjs.from_fitsfile(std_outfile)
    this_det = sobjs.DET == det
    # make the get_std from pypeit a utility function or class method
    if np.any(this_det):
        sobjs_det = sobjs[this_det]
        sobjs_std = sobjs_det.get_std()
        std_trace = None if sobjs_std is None else sobjs_std.TRACE_SPAT.flatten()
    else:
        std_trace=None

    # Build Science image
    sciImg = buildimage.buildimage_fromlist(spectrograph, det, parset['scienceframe'], sci_files, bpm=msbpm,
                                            slits=slits, ignore_saturation=False)

    # Background Image?
    sciImg = sciImg.sub(buildimage.buildimage_fromlist(spectrograph, det,  parset['scienceframe'], bg_files, bpm=msbpm,
                                                       slits=slits,ignore_saturation=False), parset['scienceframe']['process'])
    # Build the Calibrate object
    caliBrate = calibrations.Calibrations(None, parset['calibrations'], spectrograph, None)
    caliBrate.slits = slits
    caliBrate.wavetilts = tilts_obj
    caliBrate.wv_calib = wv_calib

    # Instantiate Reduce object
    # Required for pypeline specific object
    # At instantiaton, the fullmask in self.sciImg is modified
    redux = reduce.Reduce.get_instance(sciImg, spectrograph, parset, caliBrate, 'science', ir_redux=True, show=pargs.show,
                                       det=det, std_outfile=std_outfile)

    manual_extract_dict = None
    skymodel, objmodel, ivarmodel, outmask, sobjs, waveImg, tilts = redux.run(
        std_trace=std_trace, manual_extract_dict=manual_extract_dict, show_peaks=pargs.show)

    # TODO -- Do this upstream
    # Tack on detector
    for sobj in sobjs:
        sobj.DETECTOR = sciImg.detector


    # Construct the Spec2DObj
    spec2DObj = spec2dobj.Spec2DObj(det=det,
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
    spec2DObj.process_steps = sciImg.process_steps
    all_spec2d = spec2dobj.AllSpec2DObj()
    all_spec2d['meta']['ir_redux'] = True
    all_spec2d[det] = spec2DObj
    # TODO -- Should we reset/regenerate self.slits.mask for a new exposure
    save_exposure(fitstbl, 0, spectrograph, science_path, parset, caliBrate, all_spec2d, sobjs)
    embed()


    #
    # TODO -- Get the type_bits from  'science'
    bm = framematch.FrameTypeBitMask()
    file_bits = np.zeros(2, dtype=bm.minimum_dtype())
    file_bits[0] = bm.turn_on(file_bits[0], ['arc', 'science', 'tilt'])
    file_bits[1] = bm.turn_on(file_bits[0], ['arc', 'science', 'tilt'])

    ps.fitstbl.set_frame_types(file_bits)
    ps.fitstbl.set_combination_groups()
    # Extras
    ps.fitstbl['setup'] = 'A'
    # A-B
    ps.fitstbl['bkg_id'] = [2,1]


    # Config the run
    cfg_lines = ['[rdx]']
    cfg_lines += ['    spectrograph = {0}'.format('keck_mosfire')]
    cfg_lines += ['    redux_path = {0}'.format(os.path.join(os.getcwd(),'keck_mosfire_A'))]
    # Calibrations
    cfg_lines += ['[baseprocess]']
    cfg_lines += ['    use_pixelflat = False']
    cfg_lines += ['    use_illumflat = False']
    cfg_lines += ['[calibrations]']
    cfg_lines += ['    master_dir = {0}'.format(master_dir)]
    cfg_lines += ['    raise_chk_error = False']
    cfg_lines += ['[scienceframe]']
    cfg_lines += ['    [[process]]']
    cfg_lines += ['        mask_cr = False']
    cfg_lines += ['[reduce]']
    cfg_lines += ['    [[extraction]]']
    cfg_lines += ['        skip_optimal = True']
    if pargs.box_radius is not None: # Boxcar radius
        cfg_lines += ['        boxcar_radius = {0}'.format(pargs.box_radius)]
    cfg_lines += ['    [[findobj]]']
    cfg_lines += ['        skip_second_find = True']

    # Write
    ofiles = ps.fitstbl.write_pypeit('', configs=['A'], write_bkg_pairs=True, cfg_lines=cfg_lines)
    if len(ofiles) > 1:
        msgs.error("Bad things happened..")

    # Instantiate the main pipeline reduction object
    pypeIt = pypeit.PypeIt(ofiles[0], verbosity=2,
                           reuse_masters=True, overwrite=True,
                           logname='nires_proc_AB.log', show=False)
    # Run
    pypeIt.reduce_all()
    msgs.info('Data reduction complete')
    # QA HTML
    msgs.info('Generating QA HTML')
    pypeIt.build_qa()

    return 0

