"""
Script for fluxing PYPEIT 1d spectra

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from IPython import embed

from astropy.io import fits

from pypeit import log
from pypeit import PypeItError
from pypeit import inputfiles
from pypeit.spectrographs.util import load_spectrograph
from pypeit import fluxcalibrate
from pypeit.par import pypeitpar
from pypeit.scripts import scriptbase
from pypeit.sensfilearchive import SensFileArchive


class FluxCalib(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Flux calibrate 1D spectra produced by PypeIt',
                                    width=width, formatter=scriptbase.SmartFormatter,
                                    default_log_file=True)

        parser.add_argument("flux_file", type=str,
                            help="R|File to guide fluxing process.  This file must have the "
                                 "following format: \n\n"
                                 "F|flux read\n"
                                 "F|     filename | sensfile\n"
                                 "F|  spec1dfile1 | sensfile1\n"
                                 "F|  spec1dfile2 | \n"
                                 "F|     ...    \n"
                                 "F|flux end\n"
                                 "\nOR\n\n"
                                 "F|flux read\n"
                                 "F|     filename | sensfile\n"
                                 "F|  spec1dfile1 | sensfile1\n"
                                 "F|  spec1dfile2 | sensfile2\n"
                                 "F|  spec1dfile3 | sensfile3\n"
                                 "F|     ...    \n"
                                 "F|flux end\n"
                                 "\nOR\n\n"
                                 "F|[fluxcalib]\n"
                                 "F|  use_archived_sens = True\n"
                                 "F|flux read\n"
                                 "F|     filename\n"
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
        parser.add_argument("--par_outfile", default='fluxing.par', action="store_true",
                            help="Output to save the parameters")
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @classmethod
    def main(cls, args):
        """ Runs fluxing steps
        """
        # Initialize the log
        cls.init_log(args)

        # Set whether or not to check datamodel versions
        chk_version = not args.try_old

        # Load the file
        fluxFile = inputfiles.FluxFile.from_file(args.flux_file)

        # Read in spectrograph from spec1dfile header
        header = fits.getheader(fluxFile.filenames[0])
        spectrograph = load_spectrograph(header['PYP_SPEC'], pypeit_fits=True)

        # Parameters
        spectrograph_def_par = spectrograph.default_pypeit_par()
        par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                                 merge_with=(fluxFile.cfg_lines,))
        # Write the par to disk
        print("Writing the parameters to {}".format(args.par_outfile))
        par.to_config(args.par_outfile)

        # Chck the sizes of the inputs
        nspec = len(fluxFile.filenames)

        # Archived solution?
        if len(fluxFile.sensfiles) > 0: 
            sensfiles = fluxFile.sensfiles
        elif len(fluxFile.sensfiles) == 0 and par['fluxcalib']['use_archived_sens'] == True:
            # No sensfile specified, but an archived sensfunc can be used.
            sf_archive = SensFileArchive.get_instance(spectrograph.name)
            sensfiles = nspec*[sf_archive.get_archived_sensfile(fluxFile.filenames[0])]
        else:
            raise PypeItError(
                'Invalid format for .flux file.\n'
                'You must specify a single sensfile on the first line of the flux block,\n'
                'or specify a  sensfile for every spec1dfile in the flux block,\n'
                'or specify "use_archived_sens = True" to use an archived sensfile.\n'
                'Run pypeit_flux_calib --help for information on the format'
            )

        # Instantiate
        fluxcalibrate.flux_calibrate(fluxFile.filenames, sensfiles, par=par['fluxcalib'],
                                     chk_version=chk_version)
        log.info('Flux calibration complete')
        return 0

