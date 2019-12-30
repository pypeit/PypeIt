#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script runs PypeIt on a set of MultiSlit images
"""
import argparse

from pypeit import msgs

import warnings

def parser(options=None):

    parser = argparse.ArgumentParser(description='Script to run PypeIt in QuickLook on a set of MOS files')
    parser.add_argument('spectrograph', type=str, help='Name of spectograph, e.g. shane_kast_blue')
    parser.add_argument('full_rawpath', type=str, help='Full path to the raw files')
    parser.add_argument('arc', type=str, help='Arc frame')
    parser.add_argument('flat', type=str, help='Flat frame')
    parser.add_argument('science', type=str, help='Science frame')
    parser.add_argument('-b', '--box_radius', type=float, help='Set the radius for the boxcar extraction (arcsec)')

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):

    import os
    import numpy as np

    from IPython import embed

    from pypeit import pypeit
    from pypeit import pypeitsetup
    from pypeit.core import framematch

    spec = pargs.spectrograph

    # Setup
    data_files = [os.path.join(pargs.full_rawpath, pargs.arc),
                  os.path.join(pargs.full_rawpath,pargs.flat),
                  os.path.join(pargs.full_rawpath,pargs.science)]
    ps = pypeitsetup.PypeItSetup(data_files, path='./', spectrograph_name=spec)
    ps.build_fitstbl()
    # TODO -- Get the type_bits from  'science'
    bm = framematch.FrameTypeBitMask()
    bits = [bm.bits[iftype] for iftype in ['arc', 'pixelflat', 'trace', 'science', 'tilt']]
    ps.fitstbl.set_frame_types(np.array([2**bits[0]+2**bits[4], 2**bits[1] + 2**bits[2], 2**bits[3]]))  # 1=arc, 16=pixelflat, 32=science, trace=128
    ps.fitstbl.set_combination_groups()
    # Extras
    ps.fitstbl['setup'] = 'A'

    # Config the run
    cfg_lines = ['[rdx]']
    cfg_lines += ['    spectrograph = {0}'.format(spec)]
    cfg_lines += ['    redux_path = {0}_A'.format(os.path.join(os.getcwd(),spec))]
    cfg_lines += ['[calibrations]']
    cfg_lines += ['    [[scienceframe]]']
    cfg_lines += ['        [[process]]']
    cfg_lines += ['              cr_reject = False']
    cfg_lines += ['[scienceimage]']
    cfg_lines += ['    [[extraction]]']
    cfg_lines += ['         boxcar_only = True']
    if pargs.box_radius is not None: # Boxcar radius
        cfg_lines += ['    boxcar_radius = {0}'.format(pargs.box_radius)]
    cfg_lines += ['    [[findobj]]']
    cfg_lines += ['         skip_second_find = True']

    # Write
    ofiles = ps.fitstbl.write_pypeit('', configs=['A'], write_bkg_pairs=True, cfg_lines=cfg_lines)
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

