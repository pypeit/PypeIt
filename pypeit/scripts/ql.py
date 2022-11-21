"""
Script for quick-look reductions for Multislit observations.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import time

import numpy as np

from astropy.table import Table

from pypeit import utils
from pypeit.scripts import run_pypeit
from pypeit import par, msgs
from pypeit import pypeitsetup
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core import setup
from pypeit.core import quicklook
from pypeit.scripts import scriptbase
from pypeit.spectrographs import available_spectrographs

from IPython import embed

class QL(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Script to produce quick-look multislit PypeIt reductions', width=width)
        parser.add_argument('spectrograph', type=str,
                            help='A valid spectrograph identifier: {0}'.format(
                                 ', '.join(available_spectrographs)))
        parser.add_argument('--rawfile_list', type=str, 
                            help='File providing raw files to reduce including their path(s)')
        parser.add_argument('--full_rawpath', type=str, 
                            help='Full path to the raw files. Used with --rawfiles or --extension')
        parser.add_argument('--raw_extension', type=str, default='.fits',
                            help='Extension for raw files in full_rawpath.  Only use if --rawfile_list and --rawfiles are not provided')
        parser.add_argument('--rawfiles', type=str, nargs='+',
                            help='space separated list of raw frames e.g. img1.fits img2.fits.  These must exist within --full_rawpath')
        parser.add_argument('--configs', type=str, default='A',
                            help='Configurations to reduce [A,all]')
        parser.add_argument('--sci_files', type=str, 
                            help='comma separated list of raw frames to be specified as science exposures (over-rides PypeIt frame typing)')
        parser.add_argument('--spec_samp_fact', default=1.0, type=float,
                            help='Make the wavelength grid finer (spec_samp_fact < 1.0) or '
                                 'coarser (spec_samp_fact > 1.0) by this sampling factor, i.e. '
                                 'units of spec_samp_fact are pixels.')
        parser.add_argument('--spat_samp_fact', default=1.0, type=float,
                            help='Make the spatial grid finer (spat_samp_fact < 1.0) or coarser '
                                 '(spat_samp_fact > 1.0) by this sampling factor, i.e. units of '
                                 'spat_samp_fact are pixels.')
        parser.add_argument("--bkg_redux", default=False, action='store_true',
                            help='If set the script will perform difference imaging quicklook. Namely it will identify '
                                 'sequences of AB pairs based on the dither pattern and perform difference imaging sky '
                                 'subtraction and fit for residuals')
        parser.add_argument("--flux", default=False, action='store_true',
                            help='This option will multiply in sensitivity function to obtain a '
                                 'flux calibrated 2d spectrum')
        parser.add_argument("--mask_cr", default=False, action='store_true',
                            help='This option turns on cosmic ray rejection. This improves the '
                                 'reduction but doubles runtime.')
        parser.add_argument("--writefits", default=False, action='store_true',
                            help="Write the ouputs to a fits file")
        parser.add_argument('--no_gui', default=False, action='store_true',
                            help="Do not display the results in a GUI")
        parser.add_argument('--box_radius', type=float,
                            help='Set the radius for the boxcar extraction')
        parser.add_argument('--offset', type=float, default=None,
                            help='Override the automatic offsets determined from the headers. '
                                 'Offset is in pixels.  This option is useful if a standard '
                                 'dither pattern was not executed.  The offset convention is '
                                 'such that a negative offset will move the (negative) B image '
                                 'to the left.')
        parser.add_argument("--redux_path", type=str, default=os.getcwd(),
                            help="Location where reduction outputs should be stored.")
        parser.add_argument("--calib_dir", type=str, 
                            help="Location folders of calibration reductions")
        parser.add_argument("--master_dir", type=str, 
                            help="Location of PypeIt Master files used for the reduction.")
        parser.add_argument('--maskID', type=int,
                            help='Reduce this slit as specified by the maskID value')
        parser.add_argument('--embed', default=False, action='store_true',
                            help='Upon completion embed in ipython shell')
        parser.add_argument("--show", default=False, action="store_true",
                            help='Show the reduction steps. Equivalent to the -s option when '
                                 'running pypeit.')
        parser.add_argument('--det', type=str, help='Detector(s) to reduce.')
        return parser


    @staticmethod
    def main(args):

        tstart = time.time()

        # Load up the spectrograph
        spectrograph = load_spectrograph(args.spectrograph)

        # Ingest Files 
        files = setup.grab_rawfiles(
            raw_paths=[args.full_rawpath], 
            file_of_files=args.rawfile_list, 
            list_of_files=args.rawfiles) 

        # Run PypeIt Setup
        ps = pypeitsetup.PypeItSetup.from_rawfiles(files,
                                        args.spectrograph) 
        ps.run(setup_only=True, no_write_sorted=True)

        '''
        # Read in the spectrograph, config the parset
        spectrograph = load_spectrograph(args.spectrograph)
        spectrograph_cfg_lines = spectrograph.config_specific_par(files[0]).to_config()
        parset = par.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines,
                                              merge_with=config_lines(args))
        _det = parse_det(args.det, spectrograph)

        target = spectrograph.get_meta_value(files[0], 'target')
        mjds = np.zeros(nfiles)
        for ifile, file in enumerate(files):
            mjds[ifile] = spectrograph.get_meta_value(file, 'mjd', ignore_bad_header=True,
                                                      no_fussing=True)
        files = files[np.argsort(mjds)]
        '''

        # Generate PypeIt files (and folders)
        # Calibs
        if args.master_dir is None:
            calib_dir = args.calib_dir if args.calib_dir is not None else args.redux_path
            calib_pypeit_files = quicklook.generate_calib_pypeit_files(
                ps, calib_dir,
                det=args.det, configs=args.configs)

            # Process them
            quicklook.process_calibs(calib_pypeit_files)

        # Science files                                
        if args.sci_files is not None:
            sci_files = args.sci_files.split(',')
            # WORK ON THIS
            embed(header='434 of ql multi')
        else:
            sci_idx = ps.fitstbl['frametype'] == 'science'

        if np.sum(sci_idx) == 0:
            msgs.error('No science frames found in the provided files.  Add at least one or specify using --sci_files.')

        # Loop on science files to setup PypeIt file and calibs
        ps_sci_list, sci_setups, full_scifiles = [], [], []
        for dir_path, sci_file in zip(ps.fitstbl['directory'][sci_idx],
            ps.fitstbl['filename'][sci_idx]):
            # Science file and setup
            full_scifile = os.path.join(dir_path, sci_file)
            ps_sci = pypeitsetup.PypeItSetup.from_file_root(
                full_scifile, spectrograph.name, extension='')
            ps_sci.run(setup_only=True, no_write_sorted=True)

            # Calibs
            if args.master_dir is None:
                calib_pypeit_file, sci_setup =\
                    quicklook.match_science_to_calibs(
                    full_scifile, ps_sci,
                    spectrograph, calib_dir)
            else:
                print("NEED TO GRAB THE SETUP")
                embed(header='458 of ql multi')
            # Save
            ps_sci_list.append(ps_sci)
            sci_setups.append(sci_setup)
            full_scifiles.append(full_scifile)

        # Only 1 setup?
        if len(np.unique(sci_setups)) != 1:
            dtbl = Table()
            dtbl['sci_files'] = ps.fitstbl['filename'][sci_idx]
            dtbl['setup'] = sci_setups
            print(dtbl)
            msgs.error('Your science files have multiple setups.  This is not supported. Remove one more of them.')

        # Let's build the PypeIt file and link to Masters
        if args.master_dir is None:
            sci_pypeit_file, sci_pypeitFile = \
                quicklook.generate_sci_pypeitfile(
                calib_pypeit_file, 
                args.redux_path,
                full_scifiles, ps_sci_list,
                maskID=args.maskID)
        else:
            print("NEED TO GENERATE FROM SCRATCH")
            embed(header='479 of ql multi')
        
        # Run it
        redux_path = os.path.dirname(sci_pypeit_file)  # Path to PypeIt file
        run_pargs = run_pypeit.RunPypeIt.parse_args(
            [sci_pypeit_file, '-r={}'.format(redux_path)])
        run_pypeit.RunPypeIt.main(run_pargs)
        msgs.info(utils.get_time_string(time.time()-tstart))
