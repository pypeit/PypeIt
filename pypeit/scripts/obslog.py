"""
This script generates an "observing log" for a directory with a set of files
to be reduced by PypeIt.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import time
import os

from IPython import embed

import numpy as np

from pypeit.spectrographs import available_spectrographs
from pypeit.scripts import scriptbase


class ObsLog(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Construct an observing log for a set of files '
                                                 'from the provided spectrograph using '
                                                 'PypeItMetaData.', width=width)
        parser.add_argument('spec', type=str,
                            help='A valid spectrograph identifier: {0}'.format(
                                        ', '.join(available_spectrographs)))
        parser.add_argument('-r', '--root', default='current working directory', type=str,
                            help='Root to search for data files.  You can provide the top-level '
                                 'directory  (e.g., /data/Kast) or the search string up through '
                                 'the wildcard (.e.g, /data/Kast/b).  Use the --extension option '
                                 'to set the types of files to search for.')
        parser.add_argument('-k', '--keys', default=False, action='store_true',
                            help='Do not produce the log; simply list the pypeit-specific '
                                 'metadata keys available for this spectrograph and their '
                                 'associated header cards.  Metadata keys with header cards that '
                                 'are None have no simple mapping between keyword and header '
                                 'card.')
        parser.add_argument('-c', '--columns', default='pypeit', type=str,
                            help='A comma-separated list of columns to include in the output '
                                 'table.  Each column must be a valid pypeit metadata keyword '
                                 'specific to this spectrograph (run pypeit_obslog with the -k '
                                 'argument to see the valid list).  Additional valid keywords '
                                 'are directory, filename, frametype, framebit, setup, calib, '
                                 'and calibbit. If \'all\', all columns collected for the pypeit '
                                 'metadata table are included.  If \'pypeit\', the columns are '
                                 'the same as those included in the pypeit file.')
        parser.add_argument('-b', '--bad_frames', default=False, action='store_true',
                            help='Clean the output of bad frames that cannot be reduced by pypeit.')
        parser.add_argument('-t', '--bad_types', default='keep', type=str,
                            help='Dictates how frames that could not be given a valid type should '
                                 'be treated.  Options are: "keep" to include them in the output, '
                                 '"rm" to remove them from the output, "only" to only include the '
                                 'frames with unknown types in the output (i.e, the frames with '
                                 'determined types are excluded).')
        parser.add_argument('-g', '--groupings', default=True, action='store_false',
                            help='Use this option to only determine the frame type.  By default, '
                                 'the script groups frames into expected configuration and '
                                 'calibration groups, and it adds the default combination groups.')
        parser.add_argument('-i', '--interact', default=False, action='store_true',
                            help='Once the metadata table is created, start an embedded IPython '
                                 'session that you can use to interact with the table (an '
                                 'Astropy.Table called fitstbl) directly.')
        parser.add_argument('-s', '--sort', default='mjd', type=str,
                            help='Metadata keyword (pypeit-specific) to use to sort the output '
                                 'table.')
        parser.add_argument('-e', '--extension', default='.fits',
                            help='File extension; compression indicators (e.g. .gz) not required.')
        parser.add_argument('-d', '--output_path', default='current working directory',
                            help='Path to top-level output directory.')
        parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                            help='Overwrite any existing files/directories')
        parser.add_argument('-f', '--file', default=None, type=str,
                            help='Name for the ascii output file.  Any leading directory path is '
                                 'stripped; use -d to set the output directory.  If None, the '
                                 'table is just printed to stdout.  If set to \'default\', the '
                                 'file is set to [spectrograph].obslog.  Note the file will *not* '
                                 'be written if you also include the -i option to embed and '
                                 'interact with the table (you can write the table using the '
                                 'astropy.table.Table.write method in the embedded IPython '
                                 'session).  The table is always written in ascii format using '
                                 'format=ascii.fixed_with for the call to '
                                 'Astropy.table.Table.write .')
        return parser

    @staticmethod
    def main(args):

        from pypeit.spectrographs.util import load_spectrograph
        from pypeit.pypeitsetup import PypeItSetup

        # Check that input spectrograph is supported
        if args.spec not in available_spectrographs:
            raise ValueError('Instrument \'{0}\' unknown to PypeIt.\n'.format(args.spec)
                             + '\tOptions are: {0}\n'.format(', '.join(available_spectrographs))
                             + '\tSelect an available instrument or consult the documentation '
                             + 'on how to add a new instrument.')

        if args.keys:
            # Only print the metadata to header card mapping
            load_spectrograph(args.spec).meta_key_map()
            return

        if args.bad_types not in ['keep', 'rm', 'only']:
            raise ValueError(f'{args.bad_types} is not a valid keyword for the --bad_types '
                             f'argument.')

        # Generate the metadata table
        ps = PypeItSetup.from_file_root(args.root, args.spec, 
                                        extension=args.extension)
        ps.run(setup_only=True,  # This allows for bad headers
               groupings=args.groupings,
               clean_config=args.bad_frames)

        # Check the file can be written (this is here because the spectrograph
        # needs to be defined first)
        _file = args.file
        if _file == 'default':
            _file = f'{ps.spectrograph.name}.obslog'
        if _file is not None:
            _odir, _file = os.path.split(_file)
            _file = os.path.join(args.output_path, _file)
            if not os.path.isdir(args.output_path):
                os.makedirs(args.output_path)
            if not args.interact and os.path.isfile(_file) and not args.overwrite:
                raise FileExistsError(f'{_file} already exists.  Use -o to overwrite.')

        # Write/Print the data
        #'{0}'.format(time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())),
        header = ['Auto-generated PypeIt Observing Log',
                  '{0}'.format(time.strftime("%Y-%m-%d",time.localtime())),
                  f'Root file string: {args.root}']
        if args.bad_types == 'keep':
            nrows = len(ps.fitstbl)
            indx = np.ones(nrows, dtype=bool)
        elif args.bad_types == 'rm':
            indx = ps.fitstbl['frametype'] != 'None'
        elif args.bad_types == 'only':
            indx = ps.fitstbl['frametype'] == 'None'
        else:
            raise ValueError('CODING ERROR: Should never get here.')
        fitstbl = ps.fitstbl.write(output='table' if args.interact else _file, rows=indx,
                                   columns=args.columns, sort_col=args.sort,
                                   overwrite=args.overwrite, header=header)

        if args.interact:
            embed()


