#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script runs PypeIt on a set of MultiSlit images
"""

def parse_args(options=None, return_parser=False):
    import argparse

    parser = argparse.ArgumentParser(description='Script to run PypeIt in QuickLook on a set of '
                                                 'MOS files',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('spectrograph', type=str, help='Name of spectograph, e.g. shane_kast_blue')
    parser.add_argument('full_rawpath', type=str, help='Full path to the raw files')
    parser.add_argument('arc', type=str, help='Arc frame filename')
    parser.add_argument('flat', type=str, help='Flat frame filename')
    parser.add_argument('science', type=str, help='Science frame filename')
    parser.add_argument('-b', '--box_radius', type=float,
                        help='Set the radius for the boxcar extraction (arcsec)')
    parser.add_argument('-d', '--det', type=int, default=1,
                        help='Detector number. Cannot use with --slit_spat')
    parser.add_argument('--ignore_headers', default=False, action='store_true',
                        help='Ignore bad headers?')
    parser.add_argument('--user_pixflat', type=str,
                        help='Use a user-supplied pixel flat (e.g. keck_lris_blue)')
    parser.add_argument('--slit_spat', type=str,
                        help='Reduce only this slit on this detector DET:SPAT_ID, e.g. 1:175')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import os
    import numpy as np

    from IPython import embed

    from pypeit import pypeit
    from pypeit import pypeitsetup
    from pypeit.core import framematch
    from pypeit import msgs

    spec = args.spectrograph

    # Config the run
    cfg_lines = ['[rdx]']
    cfg_lines += ['    spectrograph = {0}'.format(spec)]
    cfg_lines += ['    redux_path = {0}_A'.format(os.path.join(os.getcwd(),spec))]
    if args.slit_spat is not None:
        msgs.info("--slit_spat provided.  Ignoring --det")
    else:
        cfg_lines += ['    detnum = {0}'.format(args.det)]
    # Restrict on slit
    if args.slit_spat is not None:
        cfg_lines += ['    slitspatnum = {0}'.format(args.slit_spat)]
    # Allow for bad headers
    if args.ignore_headers:
        cfg_lines += ['    ignore_bad_headers = True']
    cfg_lines += ['[scienceframe]']
    cfg_lines += ['    [[process]]']
    cfg_lines += ['        mask_cr = False']
    # Calibrations
    cfg_lines += ['[baseprocess]']
    cfg_lines += ['    use_biasimage = False']
    cfg_lines += ['[calibrations]']
    # Input pixel flat?
    if args.user_pixflat is not None:
        cfg_lines += ['    [[flatfield]]']
        cfg_lines += ['        pixelflat_file = {0}'.format(args.user_pixflat)]
    # Reduction restrictions
    cfg_lines += ['[reduce]']
    cfg_lines += ['    [[extraction]]']
    cfg_lines += ['         skip_optimal = True']
    # Set boxcar radius
    if args.box_radius is not None:
        cfg_lines += ['    boxcar_radius = {0}'.format(args.box_radius)]
    cfg_lines += ['    [[findobj]]']
    cfg_lines += ['         skip_second_find = True']

    # Data files
    data_files = [os.path.join(args.full_rawpath, args.arc),
                  os.path.join(args.full_rawpath, args.flat),
                  os.path.join(args.full_rawpath,args.science)]

    # Setup
    ps = pypeitsetup.PypeItSetup(data_files, path='./', spectrograph_name=spec,
                                 cfg_lines=cfg_lines)
    ps.build_fitstbl()
    # TODO -- Get the type_bits from  'science'
    bm = framematch.FrameTypeBitMask()
    file_bits = np.zeros(3, dtype=bm.minimum_dtype())
    file_bits[0] = bm.turn_on(file_bits[0], ['arc', 'tilt'])
    file_bits[1] = bm.turn_on(file_bits[1], ['pixelflat', 'trace', 'illumflat']
                                             if args.user_pixflat is None 
                                             else ['trace', 'illumflat'])
    file_bits[2] = bm.turn_on(file_bits[2], 'science')

    # PypeItSetup sorts according to MJD
    #   Deal with this
    asrt = []
    for ifile in data_files:
        bfile = os.path.basename(ifile)
        idx = ps.fitstbl['filename'].data.tolist().index(bfile)
        asrt.append(idx)
    asrt = np.array(asrt)

    # Set bits
    ps.fitstbl.set_frame_types(file_bits[asrt])
    ps.fitstbl.set_combination_groups()
    # Extras
    ps.fitstbl['setup'] = 'A'

    # Write
    ofiles = ps.fitstbl.write_pypeit(configs='A', write_bkg_pairs=True, cfg_lines=cfg_lines)
    if len(ofiles) > 1:
        msgs.error("Bad things happened..")

    # Instantiate the main pipeline reduction object
    pypeIt = pypeit.PypeIt(ofiles[0], verbosity=2,
                           reuse_masters=True, overwrite=True,
                           logname='mos.log', show=False)
    # Run
    pypeIt.reduce_all()
    msgs.info('Data reduction complete')
    # QA HTML
    msgs.info('Generating QA HTML')
    pypeIt.build_qa()

    return 0

