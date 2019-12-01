
from configobj import ConfigObj
import numpy as np
from pypeit import par, msgs
from astropy.io import fits
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
import argparse
import textwrap



def read_sensfile(ifile):
    """
    Read a PypeIt sens file, akin to a standard PypeIt file

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
          Contains spec1d_files and flux_files
          Empty if no flux block is specified

    """

    # Read in the pypeit reduction file
    msgs.info('Loading the fluxcalib file')
    lines = par.util._read_pypeit_file_lines(ifile)

    # Construct config to get spectrograph
    cfg_lines = list(lines)
    #cfg = ConfigObj(cfg_lines)
    #spectrograph_name = cfg['rdx']['spectrograph']
    #spectrograph = load_spectrograph(spectrograph_name)

    # Return
    return cfg_lines



def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("spec1dfile", type=str,
                        help="spec1d file for the standard that will be used to compute sensitivity function")
    parser.add_argument("algorithm", type=str, choices=['UVIS', 'IR'],
                        help=textwrap.dedent('''\
                        Specify the algorithm for computing the sensitivity function. The options are: 
                        
                           UVIS = Should be used for data with lambda < 7000A. No detailed model of telluric absorption 
                           but corrects for atmospheric extinction. 
                           
                           IR   = Should be used for data with lambbda > 7000A. Peforms joint fit for sensitivity 
                                  function and telluric absorption using HITRAN models'''))
    parser.add_argument("-o", "--outfile", type=str,
                        help="Ouput file for sensitivity function. If not specified, the sensitivity function will "
                             "be written out to a standard filename in the current working directory, i.e. if the "
                             "standard spec1d file is named spec1d_b24-Feige66_KASTb_foo.fits the sensfunc will be "
                             "written to sens_b24-Feige66_KASTb_foo.fits")
    parser.add_argument("-s", "--sens_file", type=str, help="Configuration file to change default sensivity function parameters")
    parser.add_argument("--debug", default=False, action="store_true", help="show debug plots?")
    #parser.add_argument("--plot", default=False, action="store_true", help="Show the sensitivity function?")
    parser.add_argument("--par_outfile", default='sensfunc.par', help="Name of outut file to save the parameters used by the fit")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args



def main(args):
    """ Executes sensitivity function computation.
    """

    # If the sensfile is passed read it in
    if args.sensfile is not None:
        cfg_lines = read_sensfile(args.sens_file)
    else:
        cfg_lines = None

    # Determine the spectrograph
    header = fits.getheader(args.spec1dfile)
    spectrograph = load_spectrograph(header['PYP_SPEC'])

    # Parameters
    spectrograph_def_par = spectrograph.default_pypeit_par()
    par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                             merge_with=cfg_lines)
    # Write the par to disk
    print("Writing the parameters to {}".format(args.par_outfile))
    par.to_config(args.par_outfile)
    # TODO JFH I would like to be able to run par['telluric'].to_config('coadd2d.par')