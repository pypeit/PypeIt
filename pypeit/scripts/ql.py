"""
Script for quick-look reductions for Multislit observations.

.. include:: ../include/links.rst
"""
import os
import time
import glob

import numpy as np


from pypeit import utils
from pypeit.scripts import run_pypeit
from pypeit import msgs
from pypeit import pypeitsetup
from pypeit import io
from pypeit import pypeit
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
                            help='Full path to the raw files. Used with --rawfiles or --raw_extension')
        parser.add_argument('--raw_extension', type=str, default='.fits',
                            help='Extension for raw files in full_rawpath.  Only use if --rawfile_list and --rawfiles are not provided')
        parser.add_argument('--rawfiles', type=str, nargs='+',
                            help='space separated list of raw frames e.g. img1.fits img2.fits.  These must exist within --full_rawpath')
        parser.add_argument('--configs', type=str, default='A',
                            help='Configurations to reduce [A,all]')
        parser.add_argument('--sci_files', type=str, nargs='+',
                            help='space separated list of raw frames to be specified as science exposures (over-rides PypeIt frame typing)')
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
        parser.add_argument("--calibs_only", default=False, action="store_true",
                            help='Reduce only the calibrations?')
        return parser


    @staticmethod
    def main(args):

        tstart = time.perf_counter()

        # Ingest Files 
        files = io.grab_rawfiles(
            raw_paths=[args.full_rawpath], 
            file_of_files=args.rawfile_list, 
            list_of_files=args.rawfiles) 

        # Run PypeIt Setup on all the files
        ps = pypeitsetup.PypeItSetup.from_rawfiles(files,
                                                args.spectrograph)
        ps.run(setup_only=True, no_write_sorted=True)

        # Calibrate, if necessasry
        calib_dir = args.calib_dir if args.calib_dir is not None else args.redux_path
        if args.master_dir is None:
            # Generate PypeIt files (and folders)
            calib_pypeit_files = ps.generate_ql_calib_pypeit_files(
                calib_dir, det=args.det, configs=args.configs)
            # Process them
            for calib_pypeit_file in calib_pypeit_files: 
                # Run me via the script
                redux_path = os.path.dirname(calib_pypeit_file)  # Path to PypeIt file
                pypeIt = pypeit.PypeIt(calib_pypeit_file,
                                       redux_path=redux_path, 
                                       calib_only=True)
                calib_dict = pypeIt.calib_all()
        if args.calibs_only:
            msgs.info("Calibrations only requested.  Exiting")
            return

        # Science files                                
        if args.sci_files is not None:
            sci_idx = np.in1d(ps.fitstbl['filename'], args.sci_files)
        else:
            sci_idx = ps.fitstbl['frametype'] == 'science'

        if np.sum(sci_idx) == 0:
            msgs.error('No science frames found in the provided files.  Add at least one or specify using --sci_files.')

        # Generate science setup object
        full_scifiles = [os.path.join(dir_path, sci_file) for dir_path, sci_file in zip(
            ps.fitstbl['directory'][sci_idx], ps.fitstbl['filename'][sci_idx])]
        ps_sci = pypeitsetup.PypeItSetup.from_rawfiles(
                full_scifiles, ps.spectrograph.name)
        ps_sci.run(setup_only=True, no_write_sorted=True)
        # Limit to a single science setup
        if len(ps_sci.fitstbl.configs.keys()) > 1:
            msgs.error('Your science files come from more than one setup. Please reduce them separately.')

        # Masters dir and their setup
        if args.master_dir is None:
            # Match to calibs
            _, master_dir = quicklook.match_science_to_calibs(
                ps_sci, calib_dir)
        else:
            master_dir = args.master_dir

        # Parse for setup
        master_files = glob.glob(os.path.join(
            master_dir, 'Master*'))
        if len(master_files) == 0:
            msgs.error('No Master files found in {:s}'.format(master_dir))
        masters_setup_and_bit = os.path.basename(
            master_files[0]).split('_')[1:3]

        # Build the PypeIt file and link to Masters
        sci_pypeit_file, _ = \
                quicklook.generate_sci_pypeitfile(
                    args.redux_path, 
                    full_scifiles, 
                    master_dir, 
                    masters_setup_and_bit, 
                    ps_sci, 
                    maskID=args.maskID, 
                    det=args.det)
        
        # Run it
        redux_path = os.path.dirname(sci_pypeit_file)  # Path to PypeIt file
        run_pargs = run_pypeit.RunPypeIt.parse_args(
            [sci_pypeit_file, '-r', redux_path])
        run_pypeit.RunPypeIt.main(run_pargs)
        msgs.info(f'Quicklook completed in {utils.get_time_string(time.perf_counter()-tstart)} seconds')
