"""
This script generates files to setup a PypeIt run

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os

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
        parser.add_argument('-r', '--root', default=os.getcwd(), type=str,
                            help='Root to search for data files.  You can provide the top-level '
                                 'directory  (e.g., /data/Kast) or the search string up through '
                                 'the wildcard (.e.g, /data/Kast/b).  Use the --extension option '
                                 'to set the types of files to search for.  Default is the '
                                 'current working directory.')
        parser.add_argument('-e', '--extension', default='.fits',
                            help='File extension; compression indicators (e.g. .gz) not required.')
        parser.add_argument('-d', '--output_path', default=os.getcwd(),
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
        parser.add_argument('-v', '--verbosity', type=int, default=2,
                            help='Level of verbosity from 0 to 2.')
        return parser

    @staticmethod
    def main(args):

        from pypeit.pypeitsetup import PypeItSetup

    #    if args.root is None:
    #        raise IOError('root is a required argument.  Use the -r, --root command-line option.')
        if args.spectrograph is None:
            raise IOError('spectrograph is a required argument.  Use the -s, --spectrograph '
                          'command-line option.')

        # Check that input spectrograph is supported
        if args.spectrograph not in available_spectrographs:
            raise ValueError('Instrument \'{0}\' unknown to PypeIt.\n'.format(args.spectrograph)
                             + '\tOptions are: {0}\n'.format(', '.join(available_spectrographs))
                             + '\tSelect an available instrument or consult the documentation '
                             + 'on how to add a new instrument.')

        # Get the output directory
        sort_dir = os.path.join(args.output_path, 'setup_files')

        # Initialize PypeItSetup based on the arguments
        ps = PypeItSetup.from_file_root(args.root, args.spectrograph, extension=args.extension,
                                        output_path=sort_dir)
        # Run the setup
        ps.run(setup_only=True, sort_dir=sort_dir, write_bkg_pairs=args.background, obslog=True)

        # Use PypeItMetaData to write the complete PypeIt file
        # TODO: Set cfg_split to 'all' by default?
        if args.cfg_split is not None:
            ps.fitstbl.write_pypeit(args.output_path, cfg_lines=ps.user_cfg,
                                    write_bkg_pairs=args.background,
                                    configs=[item.strip() for item in args.cfg_split.split(',')])



