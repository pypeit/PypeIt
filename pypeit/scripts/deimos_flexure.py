#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the viewing of a raw FITS file
"""
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
    msgs.info('Loading the coadd1d file')
    lines = par.util._read_pypeit_file_lines(ifile)
    is_config = np.ones(len(lines), dtype=bool)


    # Parse the fluxing block
    spec1dfiles = []
    objids_in = []
    s, e = par.util._find_pypeit_block(lines, 'flexure')
    if s >= 0 and e < 0:
        msgs.error("Missing 'flexure end' in {0}".format(ifile))
    elif (s < 0) or (s==e):
        msgs.error("Missing flexure read or [coadd1d] block in in {0}. Check the input format for the .coadd1d file".format(ifile))
    else:
        for ctr, line in enumerate(lines[s:e]):
            prs = line.split(' ')
            spec1dfiles.append(prs[0])
            if len(prs) > 1:
                msgs.error('Invalid format for .flex file.' + msgs.newline() +
                           'You must have specify only spec1dfiles in the block ')
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

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("flex_file", type=str,
                        help="R|File to guide fluxing process.\n"
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
    parser.add_argument("--debug", default=False, action="store_true", help="show debug plots?")
    parser.add_argument("--par_outfile", default='flexure.par', action="store_true", help="Output to save the parameters")

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
    par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(), merge_with=config_lines)
    # Write the par to disk
    print("Writing the parameters to {}".format(pargs.par_outfile))
    par.to_config(pargs.par_outfile)

    # Loop to my loop
    for filename in spec1dfiles:
        mdFlex = flexure.MultiDetFlexure(spec1dfile=filename, 
                                         PYP_SPEC=spectrograph.name)

        mdFlex.init_slits()


    #filename = hdu.filename()
    #tmp = filename.split('spec1d')
    #txt = filename.replace('.fits','.txt')
    #out_dir = os.path.join(data_dir,'dmost')
    #if not os.path.isdir(out_dir):
    #    os.mkdir(out_dir)
    #slit_table_file = os.path.join(out_dir, 'dmost'+tmp[1])


        # CREATE SLIT TABLE
        msgs.info("Generating slit table")
        slits, nslits = dmost_slit_matching.create_slit_table(hdu,data_dir,txt)

        # INITIAL SKY LINE STUFF
        msgs.info("Measuring sky lines")
        slits = measure_sky_lines(slits, nslits,hdu, orig=orig)

        # FIT SURFACES
        msgs.info("Fitting the surface")
        pmodel_m, pmodel_b,pmodel_los = fit_mask_surfaces(slits)

     
        # ADD TO TABLE
        msgs.info("Table time")
        fslits = update_flexure_fit(slits,nslits, hdu, pmodel_m, pmodel_b,pmodel_los,
                                    orig=orig)

        # REFIT FOR QA PLOTS
        msgs.info("Generate QA")
        qa_flexure_plots(data_dir,nslits,slits,fslits,hdu)

        msgs.info("Write to table")
        fslits.write(slit_table_file,overwrite=True)

    # ELSE READ IT IN
    if os.path.isfile(slit_table_file):
        fslits = Table.read(slit_table_file)
    
    print("All done!!")

    return fslits

def entry_point():
    main(parse_args())

if __name__ == '__main__':
    entry_point()