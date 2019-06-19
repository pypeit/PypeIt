#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script runs PypeIt on a set of Kast images
"""
import argparse

from pypeit import msgs

import warnings

def parser(options=None):

    parser = argparse.ArgumentParser(description='Script to run PypeIt on a set of Kast files')
    parser.add_argument('camera', type=str, help='blue or red')
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
    import sys
    import numpy as np

    from IPython import embed

    from astropy.io import fits

    from pypeit import pypeit
    from pypeit import pypeitsetup

    if pargs.camers == 'blue':
        spec = 'shane_kast_blue'
    elif pargs.camers == 'red':
        spec = 'shane_kast_red'
    else:
        msgs.error("Camera needs to be blue or red")

    # Setup
    data_files = [os.path.join(pargs.full_rawpath, pargs.fileA), os.path.join(pargs.full_rawpath,pargs.fileB)]
    ps = pypeitsetup.PypeItSetup(data_files, path='./', spectrograph_name=spec)
    ps.build_fitstbl()
    # TODO -- Get the type_bits from  'science'
    ps.fitstbl.set_frame_types(np.array([32]*2))  # 1=arc, 32=science
    ps.fitstbl.set_combination_groups()
    # Extras
    ps.fitstbl['setup'] = 'A'
    # A-B
    ps.fitstbl['bkg_id'] = [2,1]

    # Calibrations
    master_dir = os.getenv('NIRES_MASTERS')
    if master_dir is None:
        msgs.error("You need to set an Environmental variable NIRES_MASTERS that points at the Master Calibs")

    # Config the run
    cfg_lines = ['[rdx]']
    cfg_lines += ['    spectrograph = {0}'.format('keck_nires')]
    cfg_lines += ['    redux_path = {0}'.format(os.path.join(os.getcwd(),'keck_nires_A'))]
    cfg_lines += ['[calibrations]']
    cfg_lines += ['    caldir = {0}'.format(master_dir)]
    # Skip CR
    cfg_lines += ['    [[scienceframe]]']
    cfg_lines += ['        [[process]]']
    cfg_lines += ['              cr_reject = False']
    cfg_lines += ['[scienceimage]']
    cfg_lines += ['    boxcar_only = True']
    cfg_lines += ['    skip_second_find = True']
    # Boxcar radius
    if pargs.box_radius is not None:
        cfg_lines += ['    boxcar_radius = {0}'.format(pargs.box_radius)]

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

