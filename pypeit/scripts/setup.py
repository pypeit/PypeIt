"""
This script generates files to setup a PypeIt run

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import argparse
from IPython import embed

from pypeit.scripts import scriptbase
from pypeit.spectrographs import available_spectrographs


class Setup(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Parse data files to construct a pypeit file in '
                                                'preparation for reduction using \'run_pypeit\'',
                                    width=width)

        # TODO: Spectrograph should be a required argument
        parser.add_argument('-s', '--spectrograph', default=None, type=str,
                            help='A valid spectrograph identifier: {0}'.format(
                                    ', '.join(available_spectrographs)))
        parser.add_argument('-r', '--root', default='current working directory', type=str, nargs='+',
                            help='Root to search for data files.  You can provide the top-level '
                                 'directory  (e.g., /data/Kast) or the search string up through '
                                 'the wildcard (.e.g, /data/Kast/b).  Use the --extension option '
                                 'to set the types of files to search for.')
        parser.add_argument('-e', '--extension', default='.fits',
                            help='File extension; compression indicators (e.g. .gz) not required.')
        parser.add_argument('-d', '--output_path', default='current working directory',
                            help='Path to top-level output directory.')
        parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                            help='Overwrite any existing files/directories')
        parser.add_argument('-c', '--cfg_split', default=None, type=str,
                            help='Generate the PypeIt files and folders by input configuration.  '
                                 'To write all unique configurations identifed, use \'all\', '
                                 'otherwise provide the list of configuration letters; e.g., '
                                 '\'A,B\' or \'B,D,E\' or \'E\'.')
        parser.add_argument('-b', '--background', default=False, action='store_true',
                            help='Include the background-pair columns for the user to edit')
        parser.add_argument('-f', '--flexure', default=False, action='store_true',
                            help='Include the manual spatial shift (flexure) column for the user to edit')
        parser.add_argument('-m', '--manual_extraction', default=False, action='store_true',
                            help='Include the manual extraction column for the user to edit')
        parser.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Verbosity level between 0 [none] and 2 [all]. Default: 1. '
                                 'Level 2 writes a log with filename setup_YYYYMMDD-HHMM.log')
        parser.add_argument('-k', '--keep_bad_frames', default=False, action='store_true',
                            help='Keep all frames, even if they are identified as having '
                                 'bad/unrecognized configurations that cannot be reduced by '
                                 'pypeit.  (This is the opposite of the --bad_frames option in '
                                 'pypeit_obslog; i.e., you have to tell pypeit_setup to keep '
                                 'these frames, whereas you have to tell pypeit_obslog to remove '
                                 'them.')
        parser.add_argument('-G', '--gui', default=False, action='store_true',
                            help='Run setup in a GUI')        

        # NOTE: These are only used to prevent updates to some of the automated
        # document building just based on changes in the pypeit version or the
        # date.  These should *not* show up when users run `pypeit_setup -h` due
        # to the use of `help=argparse.SUPPRESS`.
        parser.add_argument('--version_override', type=str, default=None, help=argparse.SUPPRESS)
        parser.add_argument('--date_override', type=str, default=None, help=argparse.SUPPRESS)

        return parser

    @staticmethod
    def main(args):

        import time
        from pathlib import Path

        from pypeit import msgs
        from pypeit.pypeitsetup import PypeItSetup
        from pypeit.calibrations import Calibrations

        # Set the verbosity, and create a logfile if verbosity == 2
        msgs.set_logfile_and_verbosity('setup', args.verbosity)

        if args.spectrograph is None:
            if args.gui is False:
                raise IOError('spectrograph is a required argument.  Use the -s, --spectrograph '
                            'command-line option.')
        else:
            # Check that input spectrograph is supported
            if args.spectrograph not in available_spectrographs:
                raise ValueError(f'Instrument "{args.spectrograph}" unknown to PypeIt.\n'
                                 f'\tOptions are: {", ".join(available_spectrographs)}\n'
                                 '\tSelect an available instrument or consult the documentation '
                                 'on how to add a new instrument.')

        if args.gui:
            from pypeit.scripts.setup_gui import SetupGUI
            # Build up arguments to the GUI
            setup_gui_argv = ["-e", args.extension]
            if args.spectrograph is not None:
                setup_gui_argv += ["-s", args.spectrograph]

                # Pass root but only if there's a spectrograph, because
                # root has a default value but can't be acted upon by the GUI
                # without a spectrograph.
                if isinstance(args.root,list):
                    root_args = args.root
                else:
                    # If the root argument is a single string, convert it to a list.
                    # This can happen when the default for --root is used
                    root_args = [args.root]
                setup_gui_argv += ["-r"] + root_args
            gui_args = SetupGUI.parse_args(setup_gui_argv)
            SetupGUI.main(gui_args)

        # Initialize PypeItSetup based on the arguments
        ps = PypeItSetup.from_file_root(args.root, args.spectrograph, extension=args.extension)
        # Run the setup
        ps.run(setup_only=True, clean_config=not args.keep_bad_frames)

        # Print selected files
        output_path = Path(args.output_path).resolve()
        if args.cfg_split is None:
            # Output directory is hard-coded to be 'setup_files'
            output_path /= 'setup_files'
            if not output_path.exists():
                output_path.mkdir(parents=True)
            # Write the sorted file,
            sorted_file = output_path / ps.pypeit_file.replace('.pypeit', '.sorted')
            ps.fitstbl.write_sorted(sorted_file, write_bkg_pairs=args.background,
                                    write_manual=args.manual_extraction)
            # the calib file,
            calib_file = sorted_file.with_suffix('.calib')
            caldir = calib_file.parent / ps.par['calibrations']['calib_dir']
            Calibrations.association_summary(calib_file, ps.fitstbl, ps.spectrograph, caldir,
                                             overwrite=True)
            # and the obslog file
            obslog_file = sorted_file.with_suffix('.obslog')
            header = ['Auto-generated PypeIt Observing Log',
                      f'{0}'.format(time.strftime("%a %d %b %Y %H:%M:%S", time.localtime()))]
            ps.fitstbl.write(output=obslog_file, columns='pypeit', sort_col='mjd', overwrite=True,
                             header=header)
        else:
            if not output_path.exists():
                output_path.mkdir(parents=True)
            # Write the pypeit files.
            configs = [item.strip() for item in args.cfg_split.split(',')]
            pypeit_files = ps.fitstbl.write_pypeit(output_path=output_path, cfg_lines=ps.user_cfg, 
                                                   write_bkg_pairs=args.background,
                                                   write_manual=args.manual_extraction,
                                                   write_shift = args.flexure,
                                                   configs=configs,
                                                   version_override=args.version_override,
                                                   date_override=args.date_override)

            # Write the calib file for each written pypeit file.
            setups = [Path(p).resolve().name.split('.')[0].split('_')[-1] for p in pypeit_files]
            for i, setup in enumerate(setups):
                indx = ps.fitstbl.find_configuration(setup)
                calib_file = Path(pypeit_files[i]).resolve().with_suffix('.calib')
                caldir = calib_file.parent / ps.par['calibrations']['calib_dir']
                Calibrations.association_summary(calib_file, ps.fitstbl, ps.spectrograph,
                                                 caldir, subset=indx, overwrite=True)

