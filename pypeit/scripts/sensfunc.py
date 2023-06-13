"""
Script to determine the sensitivity function for a PypeIt 1D spectrum.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from IPython import embed

from pypeit.scripts import scriptbase


class SensFunc(scriptbase.ScriptBase):

    # TODO: Need an option here for multi_spec_det detectors, passed as a list
    # of numbers in the SensFunc parset, or as --det 3 7 on the command line
    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Compute a sensitivity function', width=width,
                                    formatter=scriptbase.SmartFormatter)
        parser.add_argument("spec1dfile", type=str,
                            help='spec1d file for the standard that will be used to compute '
                                 'the sensitivity function')
        parser.add_argument("--algorithm", type=str, default=None, choices=['UVIS', 'IR'],
                            help="R|Override the default algorithm for computing the sensitivity "
                                 "function.  Note that it is not possible to set --algorithm and "
                                 "simultaneously use a .sens file with the --sens_file option. If "
                                 "you are using a .sens file, set the algorithm there via:\n\n"
                                 "F|    [sensfunc]\n"
                                 "F|         algorithm = IR\n"
                                 "\nThe algorithm options are:\n\n"
                                 "UVIS = Should be used for data with lambda < 7000A.  No "
                                 "detailed model of telluric absorption but corrects for "
                                 "atmospheric extinction.\n\n"
                                 "IR = Should be used for data with lambbda > 7000A. Performs "
                                 "joint fit for sensitivity function and telluric absorption "
                                 "using HITRAN models.\n\n")
        parser.add_argument("--multi", type=str,
                            help="R|List of detector numbers to splice together for instruments "
                                 "with multiple detectors arranged in the spectral direction, "
                                 "e.g. --multi = '3,7'.  Note that it is not possible to set "
                                 "--multi and simultaneously use a .sens file with the "
                                 "--sens_file option.  If you are using a .sens file, set the "
                                 "multi_spec_det param there via:\n\n"
                                 "F|    [sensfunc]\n"
                                 "F|        multi_spec_det = 3,7\n"
                                 "\n")
        parser.add_argument("-o", "--outfile", type=str,
                            help='Output file for sensitivity function. If not specified, the '
                                 'sensitivity function will be written out to a standard filename '
                                 'in the current working directory, i.e. if the standard spec1d '
                                 'file is named spec1d_b24-Feige66_KASTb_foo.fits the sensfunc '
                                 'will be written to sens_b24-Feige66_KASTb_foo.fits. A QA file '
                                 'will also be written as sens_spec1d_b24-Feige66_KASTb_foo_QA.pdf '
                                 'and a file showing throughput plots to '
                                 'sens_spec1d_b24-Feige66_KASTb_foo_throughput.pdf. The same '
                                 'extensions for QA and throughput will be used if outfile is '
                                 'provided but with .fits trimmed off if it is in the filename.')
        parser.add_argument("-s", "--sens_file", type=str,
                            help='Configuration file with sensitivity function parameters')
        parser.add_argument("--debug", default=False, action="store_true",
                            help="show debug plots?")
        parser.add_argument("--par_outfile", default='sensfunc.par',
                            help="Name of output file to save the parameters used by the fit")
        parser.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Verbosity level between 0 [none] and 2 [all]. Default: 1. '
                                 'Level 2 writes a log with filename sensfunc_YYYYMMDD-HHMM.log')
        return parser

    @staticmethod
    def main(args):
        """Executes sensitivity function computation."""

        import os

        from astropy.io import fits

        from pypeit import msgs
        from pypeit import inputfiles
        from pypeit import io
        from pypeit.par import pypeitpar
        from pypeit import sensfunc
        from pypeit.spectrographs.util import load_spectrograph

        # Set the verbosity, and create a logfile if verbosity == 2
        msgs.set_logfile_and_verbosity('sensfunc', args.verbosity)

        # Check parameter inputs
        if args.algorithm is not None and args.sens_file is not None:
            msgs.error("It is not possible to set --algorithm and simultaneously use a .sens "
                       "file via the --sens_file option. If you are using a .sens file set the "
                       "algorithm there via:\n"
                       "\n"
                       "    [sensfunc]\n"
                       "         algorithm = IR\n"
                       "\n")

        if args.multi is not None and args.sens_file is not None:
            msgs.error("It is not possible to set --multi and simultaneously use a .sens file via "
                       "the --sens_file option. If you are using a .sens file set the detectors "
                       "there via:\n"
                       "\n"
                       "         [sensfunc]\n"
                       "              multi_spec_det = 3,7\n"
                       "\n")

        # Determine the spectrograph
        hdul = io.fits_open(args.spec1dfile)
        spectrograph = load_spectrograph(hdul[0].header['PYP_SPEC'])
        spectrograph_config_par = spectrograph.config_specific_par(hdul)

        # Construct a primary FITS header that includes the spectrograph's
        #   config keys for inclusion in the output sensfunc file
        primary_hdr = io.initialize_header()
        add_keys = (
            ['PYP_SPEC', 'DATE-OBS', 'TELESCOP', 'INSTRUME', 'DETECTOR']
            + spectrograph.configuration_keys() + spectrograph.raw_header_cards()
        )
        for key in add_keys:
            if key.upper() in hdul[0].header.keys():
                primary_hdr[key.upper()] = hdul[0].header[key.upper()]

        # If the .sens file was passed in read it and overwrite default parameters

        if args.sens_file is not None:
            sensFile = inputfiles.SensFile.from_file(args.sens_file)
            # Read sens file
            par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_config_par.to_config(),
                merge_with=(sensFile.cfg_lines,))
        else:
            par = spectrograph_config_par 

        # If algorithm was provided override defaults. Note this does undo .sens
        # file since they cannot both be passed
        if args.algorithm is not None:
            par['sensfunc']['algorithm'] = args.algorithm

        # If multi was set override defaults. Note this does undo .sens file
        # since they cannot both be passed
        if args.multi is not None:
            # parse
            multi_spec_det  = [int(item) for item in args.multi.split(',')]
            par['sensfunc']['multi_spec_det'] = multi_spec_det

        # TODO Add parsing of detectors here. If detectors passed from the
        # command line, overwrite the parset values read in from the .sens file

        # Write the par to disk
        msgs.info(f'Writing the parameters to {args.par_outfile}')
        par['sensfunc'].to_config(args.par_outfile, section_name='sensfunc', include_descr=False)

        # TODO JFH I would like to be able to run only
        # par['sensfunc'].to_config('sensfunc.par') but this crashes.
        # TODO: KBW - You can do that if you override the
        # pypeit.par.parset.ParSet.to_config method in the
        # pypeit.par.pypeitpar.SensFuncPar class.

        # Parse the output filename
        outfile = (os.path.basename(args.spec1dfile)).replace('spec1d','sens') \
                        if args.outfile is None else args.outfile
        # Instantiate the relevant class for the requested algorithm
        sensobj = sensfunc.SensFunc.get_instance(args.spec1dfile, outfile, par['sensfunc'],
                                                 debug=args.debug)
        # Generate the sensfunc
        sensobj.run()
        # Write it out to a file, including the new primary FITS header
        sensobj.to_file(outfile, primary_hdr=primary_hdr, overwrite=True)

        #TODO JFH Add a show_sensfunc option here and to the sensfunc classes.


