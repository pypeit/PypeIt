#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-

from pypeit import par, msgs
from astropy.io import fits
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
import argparse
from pypeit.core import telluric
import os
from pkg_resources import resource_filename


# A trick from stackoverflow to allow multi-line output in the help:
#https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

def read_tellfile(ifile):
    """
    Read a PypeIt telluric file, akin to a standard PypeIt file

    The top is a config block that sets ParSet parameters
      The spectrograph is not required

    Args:
        ifile: str
          Name of the flux file

    Returns:
        spectrograph: Spectrograph
        cfg_lines: list
          Config lines to modify ParSet values
        flux_dict: dict
          Contains spec1d_files
    """

    # Read in the pypeit reduction file
    msgs.info('Loading the telluric file')
    lines = par.util._read_pypeit_file_lines(ifile)

    # Construct config to get spectrograph
    cfg_lines = list(lines)

    # Return
    return cfg_lines

def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse', formatter_class=SmartFormatter)
    parser.add_argument("spec1dfile", type=str,
                        help="spec1d file for the standard that will be used to compute sensitivity function")
    parser.add_argument("--algorithm", type=str, default=None, choices=['qso', 'poly'],
                        help="telluric model algorithm")
    parser.add_argument("-g", "--tell_grid", type=str, help="Telluric model grid.")
    parser.add_argument("-p", "--pca_file", type=str, help="PCA pickle file")
    parser.add_argument("-t", "--tell_file", type=str, help="Configuration file to change default telluric parameters")
    parser.add_argument("-r", "--redshift", type=float, default=None, help="Object redshift")
    parser.add_argument("-n", "--norder", type=int, default=None, help="Polynomial order")
    parser.add_argument("--debug", default=False, action="store_true", help="show debug plots?")
    parser.add_argument("--plot", default=False, action="store_true", help="Show the telluric corrected spectrum")
    parser.add_argument("--par_outfile", default='telluric.par', help="Name of outut file to save the parameters used by the fit")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):
    """ Executes telluric correction.
    """

    # Determine the spectrograph
    header = fits.getheader(args.spec1dfile)
    spectrograph = load_spectrograph(header['PYP_SPEC'])
    spectrograph_def_par = spectrograph.default_pypeit_par()

    # If the .tell file was passed in read it and overwrite default parameters
    if args.tell_file is not None:
        cfg_lines = read_tellfile(args.tell_file)
        par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                                 merge_with=cfg_lines)
    else:
        par = spectrograph_def_par

    # If args was provided override defaults. Note this does undo .tell file
    if args.algorithm is not None:
        par['tellfit']['algorithm'] = args.algorithm
    if args.pca_file is not None:
        par['tellfit']['pca_file'] = args.pca_file
    if args.redshift is not None:
        par['tellfit']['redshift'] = args.redshift
    if args.norder is not None:
        par['tellfit']['polyorder'] = args.norder

    if args.tell_grid is not None:
        par['tellfit']['tell_grid'] = args.tell_grid
    elif par['sensfunc']['IR']['telgridfile'] is not None:
        par['tellfit']['tell_grid'] = par['sensfunc']['IR']['telgridfile']
    else:
        msgs.warn('No telluric grid file given. Using {:}'.format('TelFit_MaunaKea_3100_26100_R20000.fits'))
        par['tellfit']['tell_grid'] = resource_filename('pypeit', '/data/telluric/TelFit_MaunaKea_3100_26100_R20000.fits')

    # Write the par to disk
    print("Writing the parameters to {}".format(args.par_outfile))
    par['tellfit'].to_config('telluric.par', section_name='tellfit', include_descr=False)

    # Parse the output filename
    outfile = (os.path.basename(args.spec1dfile)).replace('.fits','_tellcorr.fits')
    modelfile = (os.path.basename(args.spec1dfile)).replace('.fits','_tellmodel.fits')

    # Run the telluric fitting procedure.
    if par['tellfit']['algorithm']=='qso':
        # run telluric.qso_telluric to get the final results
        TelQSO = telluric.qso_telluric(args.spec1dfile, par['tellfit']['tell_grid'], par['tellfit']['pca_file'],
                                       par['tellfit']['redshift'], modelfile, outfile,
                                       bal_mask=par['tellfit']['bal_mask'], disp=args.plot, debug=args.debug, show=args.plot)
    elif par['tellfit']['algorithm']=='poly':
        TelPoly = telluric.poly_telluric(args.spec1dfile, par['tellfit']['tell_grid'], modelfile, outfile,
                                         polyorder=par['tellfit']['polyorder'],
                                         fit_region_min=None, fit_region_max=None, func='legendre',
                                         model='exp', mask_lyman_a=True,
                                         debug_init=args.debug, debug=args.debug, show=args.plot)

