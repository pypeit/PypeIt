"""
Script for fluxing PYPEIT 1d spectra

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit import par, msgs
from pypeit.spectrographs.util import load_spectrograph
from pypeit import fluxcalibrate
from pypeit.par import pypeitpar
from pypeit.scripts import scriptbase
from pypeit.sensfilearchive import SensFileArchive

def read_fluxfile(ifile):
    """
    Read a ``PypeIt`` flux file, akin to a standard ``PypeIt`` file.

    The top is a config block that sets ParSet parameters.  The name of the
    spectrograph is required.

    Args:
        ifile (:obj:`str`):
            Name of the flux file

    Returns:
        :obj:`tuple`:  Three objects are returned: a :obj:`list` with the
        configuration entries used to modify the relevant
        :class:`~pypeit.par.parset.ParSet` parameters, a :obj:`list` with the
        names of spec1d files to be flux-calibrated, and a :obj:`list` with the
        names of the files with the sensitivity function to use.
    """
    # Read in the pypeit reduction file
    msgs.info('Loading the fluxcalib file')
    lines = par.util._read_pypeit_file_lines(ifile)
    is_config = np.ones(len(lines), dtype=bool)


    # Parse the fluxing block
    spec1dfiles = []
    sensfiles_in = []
    s, e = par.util._find_pypeit_block(lines, 'flux')
    if s >= 0 and e < 0:
        msgs.error("Missing 'flux end' in {0}".format(ifile))
    elif (s < 0) or (s==e):
        msgs.error("Missing flux block in {0}. Check the input format for the .flux file".format(ifile))
    else:
        for ctr, line in enumerate(lines[s:e]):
            prs = line.split(' ')
            spec1dfiles.append(prs[0])
            if len(prs) > 1:
                sensfiles_in.append(prs[1])
        is_config[s-1:e+1] = False

    # Construct config to get spectrograph
    cfg_lines = list(lines[is_config])

    # Return
    return cfg_lines, spec1dfiles, sensfiles_in


class FluxCalib(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Flux calibrate 1D spectra produced by PypeIt',
                                    width=width, formatter=scriptbase.SmartFormatter)

        parser.add_argument("flux_file", type=str,
                            help="R|File to guide fluxing process.  This file must have the "
                                 "following format: \n\n"
                                 "F|flux read\n"
                                 "F|  spec1dfile1 sensfile\n"
                                 "F|  spec1dfile2\n"
                                 "F|     ...    \n"
                                 "F|flux end\n"
                                 "\nOR\n\n"
                                 "F|flux read\n"
                                 "F|  spec1dfile1 sensfile1\n"
                                 "F|  spec1dfile2 sensfile2\n"
                                 "F|  spec1dfile3 sensfile3\n"
                                 "F|     ...    \n"
                                 "F|flux end\n"
                                 "\nOR\n\n"
                                 "F|[fluxcalib]\n"
                                 "F|  use_archived_sens = True\n"
                                 "F|flux read\n"
                                 "F|  spec1dfile1\n"
                                 "F|  spec1dfile2\n"
                                 "F|  spec1dfile3\n"
                                 "F|     ...    \n"
                                 "F|flux end\n"
                                 "\n"
                                 "That is, you must specify either a sensfile for all spec1dfiles "
                                 "on the first line, specify one sensfile for each spec1dfile, or "
                                 "specify no sensfiles and use an archived one.\n"
                                 "Archived sensfiles are available for the following spectrographs: " 
                                 + ",".join(SensFileArchive.supported_spectrographs()) + "\n\n")
        parser.add_argument("--debug", default=False, action="store_true",
                            help="show debug plots?")
        parser.add_argument("--par_outfile", default='fluxing.par', action="store_true",
                            help="Output to save the parameters")
#        parser.add_argument("--plot", default=False, action="store_true",
#                            help="Show the sensitivity function?")
        return parser

    @staticmethod
    def main(args):
        """ Runs fluxing steps
        """
        # Load the file
        config_lines, spec1dfiles, sensfiles_in = read_fluxfile(args.flux_file)
        # Read in spectrograph from spec1dfile header
        header = fits.getheader(spec1dfiles[0])
        spectrograph = load_spectrograph(header['PYP_SPEC'])

        # Parameters
        spectrograph_def_par = spectrograph.default_pypeit_par()
        par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                                 merge_with=config_lines)
        # Write the par to disk
        print("Writing the parameters to {}".format(args.par_outfile))
        par.to_config(args.par_outfile)

        # Chck the sizes of the inputs
        nspec = len(spec1dfiles)
        if len(sensfiles_in) == 1:
            # Use one sens file for each spec1d
            sensfiles = nspec*sensfiles_in
        elif len(sensfiles_in) == nspec:
            # One sensfile specified / spec1d
            sensfiles = sensfiles_in
        elif len(sensfiles_in) == 0 and par['fluxcalib']['use_archived_sens'] == True:
            # No sensfile specified, but an archived sensfunc can be used.
            sf_archive = SensFileArchive.get_instance(spectrograph.name)
            sensfiles = nspec*[sf_archive.get_archived_sensfile(spec1dfiles[0])]
        else:
            msgs.error('Invalid format for .flux file.' + msgs.newline() +
                       'You must specify a single sensfile on the first line of the flux block,' + msgs.newline() +
                       'or specify a  sensfile for every spec1dfile in the flux block,' + msgs.newline() +
                       'or specify "use_archived_sens = True" to use an archived sensfile.' + msgs.newline() +
                       'Run pypeit_flux_calib --help for information on the format')

        # Instantiate
        FxCalib = fluxcalibrate.FluxCalibrate.get_instance(spec1dfiles, sensfiles,
                                                           par=par['fluxcalib'], debug=args.debug)
        msgs.info('Flux calibration complete')
        return 0

