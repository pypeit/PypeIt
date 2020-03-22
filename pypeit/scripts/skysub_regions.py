#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the user to view a 2D FITS file
and define the sky background regions interactively.
Run above the Science/ folder.
"""

import os
import argparse
import numpy as np
from configobj import ConfigObj

from astropy.table import Table

from pypeit import msgs
from pypeit import slittrace
from pypeit.core.gui.skysub_regions import SkySubGUI
from pypeit.core import procimg
from pypeit.par.util import parse_pypeit_file
from pypeit.par import PypeItPar
from pypeit.masterframe import MasterFrame
from pypeit.spectrographs.util import load_spectrograph
from pypeit.metadata import PypeItMetaData


def parser(options=None):

    parser = argparse.ArgumentParser(description='Display a Raw science image and interactively define'
                                                 'the sky regions using a GUI. Run in the same folder'
                                                 'as your .pypeit file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default=None, help='PypeIt file')
    parser.add_argument('--det', default=1, type=int, help="Detector")

    return parser.parse_args() if options is None else parser.parse_args(options)


def get_science_frame(usrdata):
    """Find all of the indices that correspond to science frames
    """
    sciidx = np.array([], dtype=np.int)
    cntr = 0
    print("\nList of science frames:")
    for tt in range(len(usrdata)):
        ftype = usrdata['frametype'][tt].split(",")
        for ff in range(len(ftype)):
            if ftype[ff] == "science":
                sciidx = np.append(sciidx, tt)
                print("  ({0:d}) {1:s}".format(cntr+1, usrdata['filename'][tt]))
                cntr += 1
                break
    # Determine which science frame the user wants
    if cntr == 1:
        msgs.info("Only one science frame listed in .pypeit file - using that frame")
        idx = sciidx[0]
    else:
        ans = ''
        while True:
            ans = input(" Which frame would you like to select [1-{0:d}]: ".format(cntr))
            try:
                ans = int(ans)
                if 1 <= ans <= cntr:
                    idx = sciidx[ans-1]
                    break
            except ValueError:
                pass
            msgs.info("That is not a valid option!")
    return idx


def main(args):

    # Load fits file
    cfg_lines, data_files, frametype, usrdata, setups = parse_pypeit_file(args.file, runtime=False)

    # Get the raw fits file name
    sciIdx = get_science_frame(usrdata)
    fname = data_files[sciIdx]

    # Load the spectrograph
    cfg = ConfigObj(cfg_lines)
    spectrograph_name = cfg['rdx']['spectrograph']
    spectrograph = load_spectrograph(spectrograph_name, ifile=data_files[sciIdx])
    msgs.info('Loaded spectrograph {0}'.format(spectrograph.spectrograph))
    spectrograph_cfg_lines = spectrograph.config_specific_par(fname).to_config()
    par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, merge_with=cfg_lines)

    # Get the master key
    file_base = os.path.basename(fname)
    ftdict = dict({file_base: 'science'})
    fitstbl = PypeItMetaData(spectrograph, par, files=[data_files[sciIdx]], usrdata=Table(usrdata[sciIdx]), strict=True)
    fitstbl.finalize_usr_build(ftdict, setups[0])
    mkey = fitstbl.master_key(0, det=args.det)

    # Load the frame data
    rawimage, _, _, datasec, _ = spectrograph.get_rawimage(fname, args.det)
    rawimage = procimg.trim_frame(rawimage, datasec < 1)
    frame = spectrograph.orient_image(rawimage, args.det)

    # Set paths
    if par['calibrations']['caldir'] == 'default':
        mdir = os.path.join(par['rdx']['redux_path'], 'Masters')
    else:
        mdir = par['calibrations']['caldir']

    if not os.path.exists(mdir):
        mdir_base = os.path.join(os.getcwd(), os.path.basename(mdir))
        msgs.warn('Master file dir: {0} does not exist. Using {1}'.format(mdir, mdir_base))
        mdir = mdir_base

    # Load the slits
    slits = slittrace.SlitTraceSet.from_master(mkey, mdir)

    # Derive an appropriate output filename
    prefix = os.path.splitext(file_base)
    if prefix[1] == ".gz":
        prefix = os.path.splitext(prefix[0])[0]
    else:
        prefix = prefix[0]
    outname = "{0:s}_skyregions.fits".format(prefix)

    # Finally, initialise the GUI
    skyreg = SkySubGUI.initialize(args.det, frame, slits, outname=outname, runtime=False,
                                  printout=True)

    # Get the results
    skyreg.get_result()
