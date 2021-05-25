#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the viewing of a raw FITS file
"""
from IPython.terminal.embed import embed
import numpy as np
import os

from pypeit import par, msgs
from pypeit.spectrographs.util import load_spectrograph
from pypeit.par import pypeitpar
from pypeit.core import flexure


def read_flexfile(ifile):
    """
    Read a PypeIt .flex file, akin to a standard PypeIt file

    The top is a config block that sets ParSet parameters

    Args:
        ifile (str):
          Name of the flux file

    Returns:
        cfg_lines (list):
          Config lines to modify ParSet values
        spec1dfiles (list):
          Contains spec1dfiles to be flexure corrected

    """
    # Read in the pypeit reduction file
    msgs.info('Loading the flexure file')
    lines = par.util._read_pypeit_file_lines(ifile)
    is_config = np.ones(len(lines), dtype=bool)

    # Parse the fluxing block
    spec1dfiles = []
    objids_in = []
    s, e = par.util._find_pypeit_block(lines, 'flexure')
    if s >= 0 and e < 0:
        msgs.error("Missing 'flexure end' in {0}".format(ifile))
    elif (s < 0) or (s == e):
        msgs.error(
            "Missing flexure read block in {0}. Check the input format for the .flex file".format(ifile))
    else:
        for ctr, line in enumerate(lines[s:e]):
            prs = line.split(' ')
            spec1dfiles.append(prs[0])
            if len(prs) > 1:
                msgs.error('Invalid format for .flex file.' + msgs.newline() +
                           'You must specify only spec1dfiles in the block ')
        is_config[s-1:e+1] = False

    # Chck the sizes of the inputs
    nspec = len(spec1dfiles)

    # Construct config to get spectrograph
    cfg_lines = list(lines[is_config])

    # Return
    return cfg_lines, spec1dfiles


def parse_args(options=None, return_parser=False):
    import argparse
    from pypeit.spectrographs import available_spectrographs

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("flex_file", type=str,
                        help="R|File to guide flexure corrections for this multi-slit mode.\n"
                             "This file must have the following format: \n"
                             "\n"
                             "flexure read\n"
                             "  spec1dfile1\n"
                             "  spec1dfile2\n"
                             "     ...    \n"
                             "     ...    \n"
                             "flexure end\n"
                             "\n"
                             "\n")
    parser.add_argument("outroot", type=str, help='Output fileroot for the flexure fits saved as FITS.')
    parser.add_argument("--clobber", default=True,
                        action="store_true", help="Clobber output files")
    parser.add_argument("--debug", default=False,
                        action="store_true", help="show debug plots?")

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):

    from astropy.io import fits

    # Load the file
    config_lines, spec1dfiles = read_flexfile(pargs.flex_file)

    # Read in spectrograph from spec1dfile header
    header = fits.getheader(spec1dfiles[0])
    spectrograph = load_spectrograph(header['PYP_SPEC'])

    # Parameters
    spectrograph_def_par = spectrograph.default_pypeit_par()
    par = pypeitpar.PypeItPar.from_cfg_lines(
        cfg_lines=spectrograph_def_par.to_config(), merge_with=config_lines)

    # Loop to my loop
    for filename in spec1dfiles:
        # Instantiate
        mdFlex = flexure.MultiSlitFlexure(s1dfile=filename)
        # Initalize 
        msgs.info("Setup")
        mdFlex.init(spectrograph, par['flexure'])

        # INITIAL SKY LINE STUFF
        msgs.info("Measuring sky lines")
        mdFlex.measure_sky_lines()

        # FIT SURFACES
        msgs.info("Fitting the surface")
        mdFlex.fit_mask_surfaces()

        # Apply
        msgs.info("Applying flexure correction")
        mdFlex.update_fit()

        # REFIT FOR QA PLOTS
        msgs.info("Generate QA")
        mask = header['TARGET'].strip()
        fnames = header['FILENAME'].split('.')
        root = mask+'_'+fnames[2]
        mdFlex.qa_plots('./', root)

        # Write
        msgs.info("Write to disk")
        mdFlex.to_file(pargs.outroot+root+'.fits',
                       overwrite=pargs.clobber)

        # Apply??

    print("All done!!")


def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
