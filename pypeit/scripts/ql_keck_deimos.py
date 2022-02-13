#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script runs PypeIt *fast* on a set of DEIMOS images

.. include:: ../include/links.rst
"""
import os
import glob
import copy
import yaml
import io
import numpy as np
import argparse

from pypeit.scripts import scriptbase

from pypeit.scripts.utils import Utilities
from pypeit.par import PypeItPar
from pypeit.par.util import parse_pypeit_file
from pypeit.pypeitsetup import PypeItSetup
from pypeit import msgs
from pypeit.scripts import run_pypeit
from pypeit.par.util import make_pypeit_file
from pypeit import utils

from IPython import embed

class QLKeckDEIMOS(scriptbase.ScriptBase):

    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'pypeit_ql_keck_deimos'

    @classmethod
    def get_parser(cls, width=None):
        import argparse
        parser = super().get_parser(description='Pypeit QL for Keck/DEIMOS',
                                    width=width, formatter=argparse.RawDescriptionHelpFormatter)
        parser = argparse.ArgumentParser(description='Script to run PypeIt in QuickLook on a set of '
                                                    'Keck/DEIMOS files',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('full_rawpath', type=str, help='Full path to the raw files')
        parser.add_argument('--redux_path', type=str, 
                            help='Path to where reduction products lie')
        parser.add_argument('--root', type=str, help='Root of filenames, eg.  DE.2018')
        parser.add_argument('--calibs_only', default=False, action='store_true',
                            help='Run on calibs only')
        parser.add_argument('--science', type=str, help='Science frame filename')
        parser.add_argument('-b', '--box_radius', type=float,
                            help='Set the radius for the boxcar extraction (arcsec)')
        parser.add_argument('-d', '--det', type=int, default=0,
                            help='Detector number. Cannot use with --slit_spat')
        parser.add_argument('--ignore_headers', default=False, action='store_true',
                            help='Ignore bad headers?')
        parser.add_argument('--user_pixflat', type=str,
                            help='Use a user-supplied pixel flat (e.g. keck_lris_blue)')
        parser.add_argument('--maskID', type=int,
                            help='Reduce this slit as specified by the maskID value')
        parser.add_argument('--slit_spat', type=str,
                            help='Reduce only this slit on this detector DET:SPAT_ID, e.g. 0:175')

        return parser

    @staticmethod
    def main(pargs):

        script_Utils = Utilities('keck_deimos')

        # Afternoon calibs
        if pargs.calibs_only:
            if pargs.root is None:
                msgs.error("Need to set --root when using --calibs_only")
            msgs.info("Processing calibrations only")
            process_calibs(pargs, script_Utils)
            return

        # Slurp the afternoon runs and grab the right PypeIt file
        calib_pypeit_file, ps_sci = get_science_setup(pargs, script_Utils)

        # Run on the science frame
        run_on_science(pargs, script_Utils, calib_pypeit_file, ps_sci)


def process_calibs(pargs, script_Utils):
    """
    Process the calibration files

    Args:
        pargs (argparse.ArgumentParser):
            Command line arguments
        script_Utils (:class:`pypeit.scripts.utilts.Utilities`):
    """
    # Slurp
    output_path = pargs.redux_path

    # Run setup
    ps, setups, indx = script_Utils.run_setup(os.path.join(pargs.full_rawpath, pargs.root),
                                              extension='.fits')
    # Restrict on detector -- May remove this
    ps.user_cfg = ['[rdx]', 'spectrograph = {}'.format(ps.spectrograph.name)]
    ps.user_cfg += ['detnum = {}'.format(pargs.det)]
    # Avoid crash in flat fielding from saturated slits
    ps.user_cfg += ['[calibrations]', '[[flatfield]]', 'saturated_slits = mask']
    # Other DEIMOS fiddling here (e.g. reduce edge_thresh value)

    # TODO -- Remove the science files!  We want calibs only

    # Write the PypeIt files
    pypeit_files = ps.fitstbl.write_pypeit(output_path=output_path,
                                          cfg_lines=ps.user_cfg,
                                          configs='all')

    # Loop on setups, rename + run calibs
    for pypeit_file, setup, ii in zip(pypeit_files, setups, indx):

        # Rename with _calib
        calib_pypeit_file = pypeit_file.replace('deimos_{}.pypeit'.format(setup),
                                                'deimos_calib_{}.pypeit'.format(setup))
        os.rename(pypeit_file, calib_pypeit_file)

        # Run me via the script
        redux_path = os.path.dirname(calib_pypeit_file)  # Path to PypeIt file
        run_pargs = run_pypeit.RunPypeIt.parse_args([calib_pypeit_file,
                                           '-r={}'.format(redux_path),
                                           '-c'])
        run_pypeit.RunPypeIt.main(run_pargs)


def get_science_setup(pargs, script_Utils):
    """
    Figure out the matching setup for the science file

    Args:
        pargs (argparse.ArgumentParser):
            Command line arguments
        script_Utils (:class:`pypeit.scripts.utilts.Utilities`):

    Returns:
        tuple: int, :class:`pypeit.pypeitsetup.PypeItSetup`

    """
    # Check file exists
    science_file = os.path.join(pargs.full_rawpath, pargs.science)
    if not os.path.isfile(science_file):
        msgs.error("Your science filename {} does not exist. Check your path".format(science_file))
    # Run setup
    ps_sci, _, indx = script_Utils.run_setup(science_file,
                                              extension='', 
                                              no_write_sorted=True)
    # Check against existing PypeIt files
    pypeit_files = glob.glob(os.path.join(pargs.redux_path, 'keck_deimos_*', 'keck_deimos_calib_*.pypeit'))
    mtch = []
    for pypeit_file in pypeit_files:
        # Read
        tmp_scriptU = Utilities('keck_deimos', pypeit_file=pypeit_file)
        # Check for a match
        match = True
        for key in tmp_scriptU.spectrograph.configuration_keys():
            if ps_sci.fitstbl[key][0] != tmp_scriptU.sdict[key]:
                match = False
        if match:
            mtch.append(pypeit_file)
    if len(mtch) != 1:
        msgs.error("Matched to more than one setup.  Inconceivable!")

    return mtch[0], ps_sci

def run_on_science(pargs, script_Utils, calib_pypeit_file, ps_sci):
    """
    Process a science frame

    Args:
        pargs (argparse.ArgumentParser):
            Command line arguments
        script_Utils:
        calib_pypeit_file (str):
        ps_sci (:class:`pypeit.pypeitsetup.PypeItSetup`):
    """
    # Load up and slurp PypeIt file
    ps = PypeItSetup.from_pypeit_file(calib_pypeit_file)
    # Parse science file info
    science_file = os.path.join(pargs.full_rawpath, pargs.science)
    science_pypeit = calib_pypeit_file.replace('calib', 'science')
    # Add to file list, it not present
    if science_file not in ps.file_list:
        ps.file_list += [science_file]
        # Add to usrdata 
        new_row = {}
        for key in ps.usrdata.keys():
            new_row[key] = ps_sci.fitstbl[key][0]
        ps.usrdata.add_row(new_row)

    # Build
    _ = ps.build_fitstbl()

    # Remove any extraneous science files in the folder
    gd_files = (ps.usrdata['filename'] == os.path.basename(science_file)) | (
        ps.usrdata['frametype'] != 'science')
    # Update
    #ps.fitstbl.table = ps.fitstbl.table[gd_files]
    #ps.file_list = np.array(ps.file_list)[gd_files].tolist()
    ps.usrdata = ps.usrdata[gd_files]

    # Generate PypeIt file
    # TODO -- Generalize this with metadata.write_pypeit()
    # Write the file
    setup_lines = yaml.dump(utils.yamlify(ps.setup_dict)).split('\n')[:-1]
    with io.StringIO() as ff:
        ps.usrdata.write(ff, format='ascii.fixed_width')
        data_lines = ff.getvalue().split('\n')[:-1]
    paths = [os.path.dirname(item) for item in ps.file_list]
    paths = np.unique(paths)

    # Add to configs
    if pargs.slit_spat is not None:
        # Remove detnum
        for kk, item in enumerate(ps.user_cfg):
            if 'detnum' in item:
                ps.user_cfg.pop(kk)
        ridx = ps.user_cfg.index('[rdx]')
        ps.user_cfg.insert(ridx+1, '    slitspatnum = {0}'.format(pargs.slit_spat))
    else:
        raise NotImplementedError('NOT READY:  118 of ql_deimos')

    # Do it
    make_pypeit_file(science_pypeit,
                     script_Utils.spectrograph.name, [],
                     cfg_lines=ps.user_cfg,
                     setup_lines=setup_lines,
                     sorted_files=data_lines, paths=paths)

    # Run me!
    redux_path = os.path.dirname(science_pypeit)  # Path to PypeIt file
    run_pargs = run_pypeit.RunPypeIt.parse_args([science_pypeit,
                                   '-r={}'.format(redux_path),
                                   ])
    run_pypeit.RunPypeIt.main(run_pargs)