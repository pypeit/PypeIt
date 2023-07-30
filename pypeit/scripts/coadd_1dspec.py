"""
Script for coadding PypeIt 1d spectra

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os

from IPython import embed

import numpy as np

from astropy.io import fits
from astropy.time import Time

from pypeit import msgs
from pypeit import inputfiles
from pypeit import coadd1d
from pypeit import inputfiles
from pypeit.par import pypeitpar
from pypeit.scripts import scriptbase
from pypeit.spectrographs.util import load_spectrograph


## TODO: This is basically the exact same code as read_fluxfile in the fluxing
## script. Consolidate them? Make this a standard method in parse or io.
#def read_coaddfile(ifile):
#    """
#    Read a ``PypeIt`` coadd1d file, akin to a standard PypeIt file.
#
#    The top is a config block that sets ParSet parameters.  The name of the
#    spectrograph is required.
#
#    Args:
#        ifile (:obj:`str`):
#            Name of the coadd file
#
#    Returns:
#        :obj:`tuple`:  Three objects are returned: a :obj:`list` with the
#        configuration entries used to modify the relevant
#        :class:`~pypeit.par.parset.ParSet` parameters, a :obj:`list` the names
#        of spec1d files to be coadded, and a :obj:`list` with the object IDs
#        aligned with each of the spec1d files.
#    """
#    # Read in the pypeit reduction file
#    msgs.info('Loading the coadd1d file')
#    lines = inputfiles.read_pypeit_file_lines(ifile)
#    is_config = np.ones(len(lines), dtype=bool)
#
#
#    # Parse the fluxing block
#    spec1dfiles = []
#    objids_in = []
#    s, e = inputfiles.InputFile.find_block(lines, 'coadd1d')
#    if s >= 0 and e < 0:
#        msgs.error("Missing 'coadd1d end' in {0}".format(ifile))
#    elif (s < 0) or (s==e):
#        msgs.error("Missing coadd1d read or [coadd1d] block in in {0}. Check the input format for the .coadd1d file".format(ifile))
#    else:
#        for ctr, line in enumerate(lines[s:e]):
#            prs = line.split(' ')
#            spec1dfiles.append(prs[0])
#            if ctr == 0 and len(prs) != 2:
#                msgs.error('Invalid format for .coadd1d file.' + msgs.newline() +
#                           'You must have specify a spec1dfile and objid on the first line of the coadd1d block')
#            if len(prs) > 1:
#                objids_in.append(prs[1])
#        is_config[s-1:e+1] = False
#
#    # Chck the sizes of the inputs
#    nspec = len(spec1dfiles)
#    if len(objids_in) == 1:
#        objids = nspec*objids_in
#    elif len(objids_in) == nspec:
#        objids = objids_in
#    else:
#        msgs.error('Invalid format for .flux file.' + msgs.newline() +
#                   'You must specify a single objid on the first line of the coadd1d block,' + msgs.newline() +
#                   'or specify am objid for every spec1dfile in the coadd1d block.' + msgs.newline() +
#                   'Run pypeit_coadd_1dspec --help for information on the format')
#    # Construct config to get spectrograph
#    cfg_lines = list(lines[is_config])
#
#    # Return
#    return cfg_lines, spec1dfiles, objids


def build_coadd_file_name(spec1dfiles, spectrograph):
    """Build the output file name for coadding.
    The filename convention is coadd1d_<target>_<instrument name>_<YYYYMMDD>.fits or
    coadd1d_<target>_<instrument name>_<YYYYMMDD>-<YYYYMMDD>.fits if the coadd included more than
    one day's worth of data. The default location of the file will be along side the first spec1d file.

    Currently instrument_name is taken from spectrograph.camera

    Returns: 
        str:  The name of the coadd output file.
    """
    mjd_list = []
    for f in spec1dfiles:
        try:
            mjd_list.append(float(fits.getheader(f)['MJD']))
        except Exception as e:
            msgs.error(f"Failed to read MJD from {f}: {e}")

    start_mjd = np.min(mjd_list)
    end_mjd = np.max(mjd_list)

    start_date_portion = Time(start_mjd, format="mjd").strftime('%Y%m%d')
    end_date_portion = Time(end_mjd, format="mjd").strftime('%Y%m%d')

    date_portion = f"{start_date_portion}_{end_date_portion}"

    instrument_name = spectrograph.camera
    target = fits.getheader(spec1dfiles[0])['TARGET']
    path = os.path.dirname(os.path.abspath(spec1dfiles[0]))
    return os.path.join(path, f'coadd1d_{target}_{instrument_name}_{date_portion}.fits')

class CoAdd1DSpec(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Coadd 1D spectra produced by PypeIt',
                                    width=width, formatter=scriptbase.SmartFormatter)

        parser.add_argument('coadd1d_file', type=str,
                            help="R|File to guide coadding process. This file must have the "
                                 "following format (see docs for further details including the use of paths): \n\n"
                                 "F|[coadd1d]\n"
                                 "F|   coaddfile='output_filename.fits' # Optional\n"
                                 "F|   sensfuncfile = 'sensfunc.fits' # Required only for Echelle\n"
                                 "\n"
                                 "F|   coadd1d read\n"
                                 "F|        filename | obj_id\n"
                                 "F|     spec1dfile1 | objid1\n"
                                 "F|     spec1dfile2 | objid2\n"
                                 "F|     spec1dfile3 | objid3\n"
                                 "F|        ...    \n"
                                 "F|   coadd1d end\n"
                                 "\n OR the coadd1d read/end block can look like \n\n"
                                 "F|  coadd1d read\n"
                                 "F|        filename | obj_id\n"
                                 "F|     spec1dfile1 | objid \n"
                                 "F|     spec1dfile2 | \n"
                                 "F|     spec1dfile3 | \n"
                                 "F|     ...    \n"
                                 "F|  coadd1d end\n"
                                 "\n"
                                 "That is the coadd1d block must be a two column list of "
                                 "spec1dfiles and objids, but you can specify only a single objid for "
                                 "all spec1dfiles on the first line\n\n"
                                 "Where: \n"
                                 "\n"
                                 "spec1dfile: full path to a PypeIt spec1dfile\n\n"
                                 "objid: the object identifier. To determine the objids inspect "
                                 "the spec1d_*.txt files or run pypeit_show_1dspec spec1dfile "
                                 "--list\n\n"
                                 "If the coaddfile is not given the output file will be placed "
                                 "in the same directory as the first spec1d file.\n\n")
        parser.add_argument("--debug", default=False, action="store_true", help="show debug plots?")
        parser.add_argument("--show", default=False, action="store_true",
                            help="show QA during coadding process")
        parser.add_argument("--par_outfile", default='coadd1d.par',
                            help="Output to save the parameters")
        parser.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Verbosity level between 0 [none] and 2 [all]. Default: 1. '
                                 'Level 2 writes a log with filename coadd_1dspec_YYYYMMDD-HHMM.log')
        #parser.add_argument("--test_spec_path", type=str, help="Path for testing")
        return parser

    @staticmethod
    def main(args):
        """ Runs the 1d coadding steps
        """
        # Set the verbosity, and create a logfile if verbosity == 2
        msgs.set_logfile_and_verbosity('coadd_1dspec', args.verbosity)

        # Load the file
        #config_lines, spec1dfiles, objids = read_coaddfile(args.coadd1d_file)
        coadd1dFile = inputfiles.Coadd1DFile.from_file(args.coadd1d_file)

        # Append path for testing
        #if args.test_spec_path is not None:
        #    spec1dfiles = [os.path.join(args.test_spec_path, ifile) for ifile in spec1dfiles]

        # Read in spectrograph from spec1dfile header
        header = fits.getheader(coadd1dFile.filenames[0])
        spectrograph = load_spectrograph(header['PYP_SPEC'])

        # Parameters
        spectrograph_def_par = spectrograph.default_pypeit_par()
        par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                                 merge_with=(coadd1dFile.cfg_lines,))
        # Write the par to disk
        print("Writing the parameters to {}".format(args.par_outfile))
        par.to_config(args.par_outfile)
        sensfile = par['coadd1d']['sensfuncfile']
        coaddfile = par['coadd1d']['coaddfile']

        # Testing?
        #if args.test_spec_path is not None:
        #    if sensfile is not None:
        #        sensfile = os.path.join(args.test_spec_path, sensfile)
        #    coaddfile = os.path.join(args.test_spec_path, coaddfile)

        if spectrograph.pypeline == 'Echelle' and sensfile is None:
            msgs.error('You must specify the sensfuncfile in the .coadd1d file for Echelle coadds')

        if coaddfile is None:
            coaddfile = build_coadd_file_name(coadd1dFile.filenames, spectrograph)

        # Instantiate
        coAdd1d = coadd1d.CoAdd1D.get_instance(coadd1dFile.filenames, 
                                               coadd1dFile.objids, 
                                               spectrograph=spectrograph,
                                               par=par['coadd1d'], sensfile=sensfile,
                                               debug=args.debug, show=args.show)
        # Run
        coAdd1d.run()
        # Save to file
        coAdd1d.save(coaddfile)
        msgs.info('Coadding complete')



