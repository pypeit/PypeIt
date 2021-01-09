"""
This script generates an "observing log" for a directory with a set of files
to be reduced by PypeIt.
"""
import time
import os
import shutil
import io       # NOTE: This is *not* pypeit.io

from IPython import embed

import numpy as np

from pypeit.spectrographs import available_spectrographs
from pypeit.spectrographs.util import load_spectrograph

def parse_args(options=None, return_parser=False):
    import argparse

    # TODO: Add argument that specifies the log file
    parser = argparse.ArgumentParser(description='Construct an observing log for a set of files '
                                                 'from the provided spectrograph using '
                                                 'PypeItMetaData.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # TODO: Make root and spectrograph required arguments
    parser.add_argument('spec', type=str,
                        help='A valid spectrograph identifier: {0}'.format(
                                    ', '.join(available_spectrographs)))
    parser.add_argument('root', type=str, help='File path+root, e.g. /data/Kast/b ')

    parser.add_argument('-k', '--keys', default=False, action='store_true',
                        help='Do not produce the log; simply list the pypeit-specific metadata '
                             'keys available for this spectrograph and their associated header '
                             'cards.  Metadata keys with header cards that are None have no '
                             'simple mapping between keyword and header card.')
    parser.add_argument('-c', '--columns', default='pypeit', type=str,
                        help='A comma-separated list of columns to include in the output table. '
                             'Each column must be a valid pypeit metadata keyword specific to '
                             'this spectrograph (run pypeit_obslog with the -k argument to see '
                             'the valid list).  Additional valid keywords are directory, '
                             'filename, frametype, framebit, setup, calib, and calibbit. '
                             'If \'all\', all columns collected for the pypeit metadata table '
                             'are included.  If \'pypeit\', the columns are the same as those '
                             'included in the pypeit file.')
    parser.add_argument('-b', '--bad_frames', default=False, action='store_true',
                        help='Clean the output of bad frames that cannot be reduced by pypeit.')
    parser.add_argument('-g', '--groupings', default=False, action='store_true',
                        help='By default, the script only attempts to determine the type of '
                             'each frame.  Use this option to further group the frames into '
                             'expected configuration and calibration groups.')
    parser.add_argument('-i', '--interact', default=False, action='store_true',
                        help='Once the metadata table is created, start an embedded IPython '
                             'session that you can use to interact with the table (an '
                             'Astropy.Table called fitstbl) directly.')
    parser.add_argument('-s', '--sort', default='mjd', type=str,
                        help='Metadata keyword (pypeit-specific) to use to sort the output table.')
    parser.add_argument('-e', '--extension', default='.fits',
                        help='File extension; compression indicators (e.g. .gz) not required.')
    parser.add_argument('-d', '--output_path', default=os.getcwd(),
                        help='Path to top-level output directory.')
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')
    parser.add_argument('-f', '--file', default=None, type=str,
                        help='Name for the output file.  Any leading directory path is stripped; '
                             'use -d to set the output directory.  If None, the table is just '
                             'printed to stdout.  If set to \'default\', the file is set to '
                             '[spectrograph].obslog.  Note the file will *not* be written if you '
                             'also include the -i option to embed and interact with the table '
                             '(you can write the table using the astropy.table.Table.write method '
                             'in the embedded IPython session).')

    # Option to write *all* or a selection of additional fitstbl columns
    # Option to list available fitstbl columns

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    from pypeit.pypeitsetup import PypeItSetup

    # Check that input spectrograph is supported
    if args.spec not in available_spectrographs:
        raise ValueError('Instrument \'{0}\' unknown to PypeIt.\n'.format(args.spec)
                         + '\tOptions are: {0}\n'.format(', '.join(available_spectrographs))
                         + '\tSelect an available instrument or consult the documentation '
                         + 'on how to add a new instrument.')

    if args.keys:
        spec = load_spectrograph(args.spec)
        meta_keys = list(spec.meta.keys())
        meta_cards = [str(spec.meta[key]['card']) for key in meta_keys]
        nk = max(12, max([len(key) for key in meta_keys]))
        nc = max(11, max([len(card) for card in meta_cards]))
        print('')
        print('{0} {1}'.format('Metadata Key'.center(nk), 'Header Card'.center(nc)))
        print('-'*nk + ' ' + '-'*nc)
        for key, card in zip(meta_keys, meta_cards):
            print('{0} {1}'.format(key.rjust(nk), card.rjust(nc)))
        print('')
        return

    # Initialize PypeItSetup based on the arguments
    # TODO: We need a way of instantiating this without requiring it to write a
    # file...
    tmpdir = os.path.join(args.output_path, 'tmplog')
    # Don't delete the directory if already existed.
    already_existed = os.path.isdir(tmpdir)
    ps = PypeItSetup.from_file_root(args.root, args.spec, extension=args.extension,
                                    output_path=tmpdir)
    if not already_existed:
        shutil.rmtree(tmpdir)
    meta_keys = list(ps.spectrograph.meta.keys())

    # Check the file can be written (this is here because the spectrograph
    # needs to be defined first)
    _file = args.file
    if _file == 'default':
        _file = '{0}.obslog'.format(ps.spectrograph.name)
    if _file is not None:
        _odir, _file = os.path.split(_file)
        if len(_odir) == 0:
            _file = os.path.join(args.output_path, _file)
        if not args.interact and os.path.isfile(_file) and not args.overwrite:
            raise FileExistsError('{0} already exists.  Use -o to overwrite.'.format(_file))

    # The following is very close to PypeitSetup.run. I've pulled out the
    # needed parts here so that I don't need to add additional keyword flags to
    # that method...

    # Initialize the metadata table
    ps.build_fitstbl(strict=False)

    # Remove frames taken instrument configurations that pypeit can't reduce
    if args.bad_frames:
        ps.fitstbl.clean_configurations()

    # Frame typing
    ps.get_frame_types(flag_unknown=True)

    if args.groupings:
        # Determine the configurations and assign each frame to the specified
        # configuration
        ps.fitstbl.set_configurations(ps.fitstbl.unique_configurations())

        # Assign frames to calibration groups
        ps.fitstbl.set_calibration_groups()

    # TODO: Much of the following should probably get stuffed into
    # PypeItMetaData.write

    # Get the columns to return
    if args.columns == 'all':
        tbl_cols = list(ps.fitstbl.keys())
    elif args.columns == 'pypeit':
        tbl_cols = ps.fitstbl.set_pypeit_cols()
    else:
        all_cols = list(ps.fitstbl.keys())
        tbl_cols = args.columns.split(',')
        badcol = [col not in all_cols for col in tbl_cols]
        if np.any(badcol):
            raise ValueError('The following columns are not valid: {0}'.format(
                             ', '.join(tbl_cols[badcol])))

    # Make sure the basic parameters are the first few columns; do them in
    # reverse order so I can always insert at the beginning of the list
    for col in ['framebit', 'frametype', 'filename', 'directory']:
        if col not in tbl_cols:
            continue
        indx = np.where([t == col for t in tbl_cols])[0][0]
        if indx != 0:
            tbl_cols.insert(0, tbl_cols.pop(indx))

    if args.sort not in ps.fitstbl.keys():
        raise ValueError(f'Cannot sort by {args.sort}.  Not a valid key.')
    ps.fitstbl.table.sort(args.sort)

    if args.interact:
        fitstbl = ps.fitstbl.table[tbl_cols]
        embed()
        return

    with io.StringIO() as ff:
        ps.fitstbl.table[tbl_cols].write(ff, format='ascii.fixed_width')
        data_lines = ff.getvalue().split('\n')[:-1]

    if _file is None:
        print('\n'.join(data_lines))
        return

    # Write the output. NOTE: Should not get here if setting overwrite=True is
    # not a valid thing to do (either the file doesn't exist or the user has
    # set the -o option)
    with open(_file, 'w') as f:
        f.write('# Auto-generated PypeIt Observing Log\n')
        f.write('# {0}\n'.format(time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
        f.write('\n')
        f.write('\n'.join(data_lines))
        f.write('\n')


