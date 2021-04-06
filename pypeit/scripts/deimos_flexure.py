#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the viewing of a raw FITS file
"""
import numpy as np

from pypeit import par, msgs

# TODO This is basically the exact same code as read_fluxfile in the fluxing script. Consolidate them? Make this
# a standard method in parse or io.
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

    parser.add_argument('file', type = str, default = None, help = 'FITS file')
    parser.add_argument('spectrograph', type=str,
                        help='A valid spectrograph identifier: {0}'.format(
                             ', '.join(available_spectrographs)))
    parser.add_argument("--list", default=False, help="List the extensions only?", action="store_true")
    parser.add_argument('--exten', type=int, default = 0, help="FITS extension")
    parser.add_argument('--det', type=int, default=1, help="Detector number (ignored for keck_lris, keck_deimos")
    parser.add_argument('--chname', type=str, default='Image', help="Name of Ginga tab")

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):

    from astropy.io import fits


    filename = hdu.filename()
    tmp = filename.split('spec1d')
    
    txt = filename.replace('.fits','.txt')
    
    out_dir = os.path.join(data_dir,'dmost')
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    slit_table_file = os.path.join(out_dir, 'dmost'+tmp[1])

    # IF FILE DOESN"T EXIST GENERATE
    if (not os.path.isfile(slit_table_file)) | (clobber == 1):

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


