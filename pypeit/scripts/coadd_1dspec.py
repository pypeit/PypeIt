#!/usr/bin/env python
"""
Script for coadding PypeIt 1d spectra
"""
from configobj import ConfigObj
import numpy as np
from pypeit import par, msgs
import argparse
from pypeit import coadd1d
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from astropy.io import fits

from IPython import embed


# A trick from stackoverflow to allow multi-line output in the help:
#https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def read_coaddfile(ifile):
    """
    Read a PypeIt .coadd1d file, akin to a standard PypeIt file

    The top is a config block that sets ParSet parameters
      The spectrograph is required

    Args:
        ifile: str
          Name of the flux file

    Returns:
        cfg_lines (list):
          Config lines to modify ParSet values
        spec1dfiles (list):
          Contains spec1dfiles to be coadded
        objids (list):
          Object ids aligned with each of the spec1dfiles


    """
    # Read in the pypeit reduction file
    msgs.info('Loading the coadd1d file')
    lines = par.util._read_pypeit_file_lines(ifile)
    is_config = np.ones(len(lines), dtype=bool)


    # Parse the fluxing block
    spec1dfiles = []
    objids = []
    sensfile = None
    s, e = par.util._find_pypeit_block(lines, 'coadd1d')
    if s >= 0 and e < 0:
        msgs.error("Missing 'coadd1d end' in {0}".format(ifile))
    elif (s < 0) or (s==e):
        msgs.error("Missing coadd1d block in in {0}. Check the input format for the .coadd1d file".format(ifile))
    else:
        for ctr, line in enumerate(lines[s:e]):
            prs = line.split(' ')
            if ctr == 0:
                # First line can have two or three entries, where third line is sensfile
                if ((len(prs) != 2) or (len(prs) != 3)):
                    msgs.error('Invalid format for .coadd1d file. First line of the coadd1d block must have format' +  msgs.newline()
                    +   'spec1dfile1 objid1' + msgs.newline()
                    +   '     OR' + msgs.newline()
                    +   'spec1dfile1 objid1 sensfile')
                elif len(prs) == 3:
                    sensfile = prs[3]
            else:
                # All other lines must have two entries
                if len(prs) != 2:
                    msgs.error('Invalid format for .coadd1d file.' + msgs.newline() +
                               'You must have specify a spec1dfile and objid on each line of the coadd1d block')
            spec1dfiles.append(prs[0])
            objids.append(prs[1])
        is_config[s-1:e+1] = False

    # Chck the sizes of the inputs
    if len(spec1dfiles) != len(objids):
        msgs.error('Problem with your .coadd1d file input')
    # Construct config to get spectrograph
    cfg_lines = list(lines[is_config])
    #cfg = ConfigObj(cfg_lines)
    #spectrograph_name = cfg['rdx']['spectrograph']
    #spectrograph = load_spectrograph(spectrograph_name)

    # Return
    return cfg_lines, spec1dfiles, objids, sensfile


def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse', formatter_class=SmartFormatter)
    parser.add_argument("coadd1d_file", type=str,
                        help="R|File to guide coadding process.\n"
                             "This file must have the following format: \n"
                             "\n"
                             "For MultiSlit:"
                             "\n"
                             "coadd1d read\n"
                             "  spec1dfile1 objid1\n"
                             "  spec1dfile2 objid2\n"
                             "  spec1dfile3 objid3\n"
                             "     ...    \n"
                             "coadd1d end\n"
                             "\n"
                             "For Echelle:"
                             "\n"
                             "  spec1dfile1 objid1 sensfile\n"
                             "  spec1dfile2 objid2\n"
                             "  spec1dfile3 objid3\n"
                             "     ...    \n"
                             "Where a spec1dfile is the path to a PypeIt spec1dfile, and objid is the object identifier.\n"
                             "To determine the objids inspect the spec1d_*.txt files or run pypeit_show_1dspec --list\n"
                             " For Echelle coadds, the file containing the sensitivity function must be passed on the first line"
                             "\n")
    parser.add_argument("--sensfile", default=None, action="store_true", help="Sensitivity function file "
                                                                              "must be passed in for Echelle coadds")
    parser.add_argument("--debug", default=False, action="store_true", help="show debug plots?")
    parser.add_argument("--show", default=False, action="store_true", help="show QA during coadding process")
    parser.add_argument("--par_outfile", default='coadd1d.par', action="store_true", help="Output to save the parameters")
#    parser.add_argument("--plot", default=False, action="store_true", help="Show the sensitivity function?")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):
    """ Runs the 1d coadding steps
    """
    # Load the file
    config_lines, spec1dfiles, objids = read_coaddfile(args.flux_file)
    # Read in spectrograph from spec1dfile header
    header = fits.getheader(spec1dfiles[0])
    spectrograph = load_spectrograph(header['PYP_SPEC'])

    # Parameters
    spectrograph_def_par = spectrograph.default_pypeit_par()
    par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(), merge_with=config_lines)
    # Write the par to disk
    print("Writing the parameters to {}".format(args.par_outfile))
    par.to_config(args.par_outfile)

    # Instantiate
    coadd = coadd1d.CoAdd1d.get_instance(spec1dfiles, sensfiles, spectrograph, par['coadd1d'], debug=args.debug, show=args.show)
    msgs.info('Coadding complete')

