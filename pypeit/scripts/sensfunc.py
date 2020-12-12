
from configobj import ConfigObj
import numpy as np
from pypeit import par, msgs
from astropy.io import fits
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
import argparse
from IPython import embed
import textwrap
from pypeit import sensfunc
import os


# A trick from stackoverflow to allow multi-line output in the help:
#https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

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



# TODO Need an option here for multi_spec_det detectors, passed as a list of numbers in the SensFunc parset, or as --det 3 7 on the command line
def parse_args(options=None, return_parser=False):
    parser = argparse.ArgumentParser(description='Compute a sensitivity function',
                                     formatter_class=SmartFormatter)
    parser.add_argument("spec1dfile", type=str,
                        help="spec1d file for the standard that will be used to compute sensitivity function")
    parser.add_argument("--algorithm", type=str, default=None, choices=['UVIS', 'IR'],
                        help="R|Override the default algorithm for computing the sensitivity function. \n"
                             "Note that it is not possible to set --algorithm and simultaneously use a .sens file with\n"
                             "the --sens_file option. If you are using a .sens file set the algorithm there via:\n"
                             "\n"
                             "    [sensfunc]\n"
                             "         algorithm = IR\n"
                             "\n"
                             "The algorithm options are:\n"
                             "\n"
                             "    UVIS = Should be used for data with lambda < 7000A.\n" 
                             "    No detailed model of telluric absorption but corrects for atmospheric extinction.\n"
                             "\n"
                             "    IR   = Should be used for data with lambbda > 7000A.\n"
                             "    Peforms joint fit for sensitivity function and telluric absorption using HITRAN models.\n"
                             "\n")
    parser.add_argument("--multi", type=str,
                        help="R|List of detector numbers to splice together for instruments with multiple detectors\n"
                             "arranged in the spectral direction, e.g. --multi = '3,7'\n"
                             "Note that it is not possible to set --multi and \n"
                             "simultaneously use a .sens file with the --sens_file option.\n"
                             "If you are using a .sens file set the multi_spec_det param there via:\n"
                             "\n"
                             "         [sensfunc]\n"
                             "              multi_spec_det = 3,7\n"
                             "\n")
    parser.add_argument("-o", "--outfile", type=str,
                        help="Ouput file for sensitivity function. If not specified, the sensitivity function will "
                             "be written out to a standard filename in the current working directory, i.e. if the "
                             "standard spec1d file is named spec1d_b24-Feige66_KASTb_foo.fits the sensfunc will be "
                             "written to sens_b24-Feige66_KASTb_foo.fits")
    parser.add_argument("-s", "--sens_file", type=str, help="Configuration file to change default sensivity function parameters")
    parser.add_argument("--debug", default=False, action="store_true", help="show debug plots?")
    parser.add_argument("--debug_init", default=False, action="store_true",
                        help="debug the initilization of the sensfunc + telluric fit for the IR algorithm")
    #parser.add_argument("--plot", default=False, action="store_true", help="Show the sensitivity function?")
    parser.add_argument("--par_outfile", default='sensfunc.par', help="Name of outut file to save the parameters used by the fit")

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
    """ Executes sensitivity function computation.
    """

    # Check parameter inputs
    if args.algorithm is not None and args.sens_file is not None:
        msgs.error("It is not possible to set --algorithm and simultaneously use a .sens file via\n"
                   "the --sens_file option. If you are using a .sens file set the algorithm there via:\n"
                   "\n"
                   "    [sensfunc]\n"
                   "         algorithm = IR\n"
                   "\n")

    if args.multi is not None and args.sens_file is not None:
        msgs.error("It is not possible to set --multi and simultaneously use a .sens file via\n"
                   "the --sens_file option. If you are using a .sens file set the detectors there via:\n"
                   "\n"
                   "\n"
                   "         [sensfunc]\n"
                   "              multi_spec_det = 3,7\n"
                   "\n")
    # Determine the spectrograph
    header = fits.getheader(args.spec1dfile)
    spectrograph = load_spectrograph(header['PYP_SPEC'])
    spectrograph_def_par = spectrograph.default_pypeit_par()
    # If the .sens file was passed in read it and overwrite default parameters
    if args.sens_file is not None:
        cfg_lines = read_sensfile(args.sens_file)
        par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                                 merge_with=cfg_lines)
    else:
        par = spectrograph_def_par

    # If algorithm was provided override defaults. Note this does undo .sens file since they cannot both be passed
    if args.algorithm is not None:
        par['sensfunc']['algorithm'] = args.algorithm
    # If multi was set override defaults. Note this does undo .sens file since they cannot both be passed
    if args.multi is not None:
        # parse
        multi_spec_det  = [int(item) for item in args.multi.split(',')]
        par['sensfunc']['multi_spec_det'] = multi_spec_det

    # TODO Add parsing of detectors here. If detectors passed from the command line, overwrite the parset values read
    # in from the .sens file

    # Write the par to disk
    print("Writing the parameters to {}".format(args.par_outfile))
    par['sensfunc'].to_config('sensfunc.par', section_name='sensfunc', include_descr=False)
    # TODO JFH I would like to be able to run only par['sensfunc'].to_config('sensfunc.par') but this crashes.

    # Parse the output filename
    outfile = (os.path.basename(args.spec1dfile)).replace('spec1d','sens') if args.outfile is None else args.outfile
    # Instantiate the relevant class for the requested algorithm
    sensobj = sensfunc.SensFunc.get_instance(args.spec1dfile, outfile, par=par['sensfunc'], debug=args.debug)
    # Generate the sensfunc
    sensobj.run()
    # Write it out to a file
    sensobj.save()

    #TODO JFH Add a show_sensfunc option here and to the sensfunc classes.