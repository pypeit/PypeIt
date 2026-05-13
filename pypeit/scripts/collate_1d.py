"""
This script collates multiple 1d spectra in multiple files by object, 
runs flux calibration on them, and then coadds them.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from datetime import datetime
import os.path
import sys

from astropy.coordinates import Angle
from astropy.time import Time

from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit import log
from pypeit.utils import is_float
from pypeit.core.collate import collate_1d
from pypeit.scripts import scriptbase
from pypeit import inputfiles



def build_parameters(args):
    """
    Read the command-line arguments and the input ``.collate1d`` file (if any), 
    to build the parameters needed by ``collate_1d``. The command line arguments
    will take priority over any parameters in the ``.collate1d`` file.

    Args:
        args (`argparse.Namespace`_):
            The parsed command line as returned by the ``argparse`` module.

    Returns:
        :obj:`tuple`: Returns three objects: a
        :class:`~pypeit.par.pypeitpar.PypeItPar` instance with the parameters
        for collate_1d, a
        :class:`~pypeit.spectrographs.spectrograph.Spectrograph` instance with
        the spectrograph parameters used to take the data, and a :obj:`list`
        with the spec1d files read from the command line or ``.collate1d`` file.
    """
    # First we need to get the list of spec1d files
    if args.input_file is not None:
        collateFile = inputfiles.Collate1DFile.from_file(args.input_file)
        cfg_lines, spec1d_files = collateFile.cfg_lines, collateFile.filenames

        # Look for a coadd1d file
        input_file_root, input_file_ext = os.path.splitext(args.input_file)
        coadd1d_config_name = input_file_root + ".coadd1d"
        if os.path.exists(coadd1d_config_name):
            coadd1DFile = inputfiles.Coadd1DFile.from_file(coadd1d_config_name)
            cfg_lines += coadd1DFile.cfg_lines

    else:
        cfg_lines = None
        spec1d_files = []

    if args.spec1d_files is not None and len(args.spec1d_files) > 0:
        spec1d_files = args.spec1d_files

    if spec1d_files is None or len(spec1d_files) == 0:
        parser = Collate1D.get_parser()
        print("Missing arguments: A list of spec1d files must be specified via command line or config file.")
        parser.print_usage()
        sys.exit(1)

    # Get the spectrograph for these files and then create a ParSet. 
    spectrograph = load_spectrograph(spec1d_files[0], pypeit_fits=True)
    spectrograph_def_par = spectrograph.default_pypeit_par()

    if cfg_lines is not None:
        # Build using config file
        params = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(), 
                                                    merge_with=(cfg_lines,))
    else:
        # No config file, use the defaults and supplement with command line args
        params = spectrograph_def_par
        params['collate1d'] = pypeitpar.Collate1DPar()

    # command line arguments take precedence over config file parameters
    if args.tolerance is not None:
        params['collate1d']['tolerance'] = args.tolerance

    if args.match_using is not None:
        params['collate1d']['match_using'] = args.match_using

    if args.exclude_slit_trace_bm is not None and len(args.exclude_slit_trace_bm) > 0:
        params['collate1d']['exclude_slit_trace_bm'] = args.exclude_slit_trace_bm.split(',')

    if args.exclude_serendip:
        params['collate1d']['exclude_serendip'] = True

    if args.wv_rms_thresh is not None:
        params['collate1d']['wv_rms_thresh'] = args.wv_rms_thresh

    if args.dry_run:
        params['collate1d']['dry_run'] = True

    if args.ignore_flux:
        params['collate1d']['ignore_flux'] = True

    if args.flux:
        params['collate1d']['flux'] = True

    if args.outdir is not None:
        params['collate1d']['outdir'] = args.outdir

    if args.spec1d_outdir is not None:
        params['collate1d']['spec1d_outdir'] = args.spec1d_outdir

    if args.refframe is not None:
        params['collate1d']['refframe'] = args.refframe

    if args.chk_version is True:
        params['rdx']['chk_version'] = True

    return params, spectrograph, spec1d_files


class Collate1D(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        # A blank Colate1DPar to avoid duplicating the help text.
        blank_pypar = pypeitpar.PypeItPar()
        blank_par = blank_pypar['collate1d']

        parser = super().get_parser(description='Flux/Coadd multiple 1d spectra from multiple '
                                                'nights and prepare a directory for the KOA.',
                                    width=width, formatter=scriptbase.SmartFormatter,
                                    default_log_file=True)

        # TODO: Is the file optional?  If so, shouldn't the first argument start
        # with '--'?
        parser.add_argument('input_file', type=str,
                            help='R|(Optional) File for guiding the collate process.  '
                                 'Parameters in this file are overidden by the command line. The '
                                 'file must have the following format:\n'
                                 '\n'
                                 'F|[collate1d]\n'
                                 'F|  tolerance             <tolerance>\n'
                                 'F|  outdir                <directory to place output files>\n'
                                 'F|  spec1d_outdir         <directory to place modified spec1ds, if any>\n'
                                 'F|  exclude_slit_trace_bm <slit types to exclude>\n'
                                 'F|  exclude_serendip      If set serendipitous objects are skipped.\n'  
                                 'F|  match_using           Whether to match using "pixel" or\n'
                                 'F|                        "ra/dec"\n'
                                 'F|  dry_run               If set the matches are displayed\n'
                                 'F|                        without any processing\n'
                                 'F|  flux                  Flux calibrate using archived sensfuncs.\n'
                                 'F|  ignore_flux           Ignore any flux calibration information in\n'
                                 'F|                        spec1d files.\n'
                                 'F|  wv_rms_thresh         If set, any objects with a wavelength rms > than the input\n'
                                 'F|                        value are skipped, else all wavelength rms values are accepted.\n'
                                 'F|  refframe              Perform reference frame correction prior to coadding.\n'
                                f'F|                        Options are {pypeitpar.WavelengthSolutionPar.valid_reference_frames()}. Defaults to None.\n'
                                 '\n'
                                 'F|spec1d read\n'
                                 'F|<path to spec1d files, wildcards allowed>\n'
                                 'F|...\n'
                                 'F|end\n',                        
                            nargs='?')
        parser.add_argument('--spec1d_files', type=str, nargs='*',
                            help='One or more spec1d files to flux/coadd/archive. '
                                 'Can contain wildcards')
        parser.add_argument('--par_outfile', default=None, type=str,
                            help='Output to save the parameters')
        parser.add_argument('--outdir', type=str, help=blank_par.descr['outdir'] + " Defaults to the current directory.")
        parser.add_argument('--spec1d_outdir', type=str, help=blank_par.descr['spec1d_outdir'] + " Defaults to overwriting existing spec1ds.")
        parser.add_argument('--tolerance', type=str, help=blank_par.descr['tolerance'])
        parser.add_argument('--match_using', type=str, choices=blank_par.options['match_using'],
                            help=blank_par.descr['match_using'])
        parser.add_argument('--dry_run', action='store_true', help=blank_par.descr['dry_run'])
        parser.add_argument('--ignore_flux', default=False, action='store_true', help=blank_par.descr['ignore_flux'])
        parser.add_argument('--flux', default=False, action = 'store_true', help=blank_par.descr['flux'])
        parser.add_argument('--exclude_slit_trace_bm', type=str, # nargs='*',
                            help=blank_par.descr['exclude_slit_trace_bm']+' Comma separated.')
        parser.add_argument('--exclude_serendip', action='store_true',
                            help=blank_par.descr['exclude_serendip'])
        parser.add_argument("--wv_rms_thresh", type=float, default = None, help=blank_par.descr['wv_rms_thresh'])
        parser.add_argument("--refframe", type=str, default = None, choices = pypeitpar.WavelengthSolutionPar.valid_reference_frames(),
                            help=blank_par.descr['refframe'])
        parser.add_argument('--chk_version', action = 'store_true', help=blank_pypar['rdx'].descr['chk_version'])
        return parser

    @classmethod
    def main(cls, args):
        # Initialize the log
        cls.init_log(args)

        start_time = datetime.now()
        (par, spectrograph, spec1d_files) = build_parameters(args)

        outdir = par['collate1d']['outdir'] 
        os.makedirs(outdir, exist_ok=True)

        # Write the par to disk
        if args.par_outfile is None:
            args.par_outfile = os.path.join(outdir, 'collate1d.par')
        print("Writing the parameters to {}".format(args.par_outfile))
        # Gather up config lines for the sections relevant to collate_1d
        config_lines = par['collate1d'].to_config(section_name='collate1d',include_descr=False) + ['']
        config_lines += par['coadd1d'].to_config(section_name='coadd1d',include_descr=False)
        if par['collate1d']['flux']:
            config_lines += [''] + par['fluxcalib'].to_config(section_name='fluxcalib',include_descr=False)
        with open(args.par_outfile, "w") as f:
            for line in config_lines:
                print (line, file=f)

        # Parse the tolerance based on the match type
        if par['collate1d']['match_using'] == 'pixel':
            tolerance = float(par['collate1d']['tolerance'])
        else:
            # For ra/dec matching, the default unit is arcseconds. We check for
            # this case by seeing if the passed in tolerance is a floating point number
            if is_float(par['collate1d']['tolerance']):
                tolerance =  float(par['collate1d']['tolerance'])
            else:
                tolerance = Angle(par['collate1d']['tolerance']).arcsec

        collate_1d(par, spectrograph, tolerance, spec1d_files)

        total_time = datetime.now() - start_time

        log.info(f'Total duration: {total_time}')

        return 0
