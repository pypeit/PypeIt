# Automated Reduction Pipeline for one arm of P200 DBSP
import os

from pypeit.pypeitsetup import PypeItSetup
from pypeit import defs
from pypeit import pypeit
from pypeit import msgs

def setup(args):
    # Check that the spectrograph is provided if using a file root
    if args['root'] is not None:
        if args['spectrograph'] is None:
            raise ValueError('Must provide spectrograph identifier with file root.')
        # Check that input spectrograph is supported
        instruments_served = defs.pypeit_spectrographs
        if args['spectrograph'] not in instruments_served:
            raise ValueError('Instrument \'{0}\' unknown to PypeIt.\n'.format(args['spectrograph'])
                             + '\tOptions are: {0}\n'.format(', '.join(instruments_served))
                             + '\tSelect an available instrument or consult the documentation '
                             + 'on how to add a new instrument.')

    # Get the output directory
    output_path = os.getcwd() if args['output_path'] is None else args['output_path']
    sort_dir = os.path.join(output_path, 'setup_files')
   
    # Initialize PypeItSetup based on the arguments
    if args['root'] is not None:
        ps = PypeItSetup.from_file_root(args['root'], args['spectrograph'], extension=args['extension'],
                                        output_path=sort_dir)
    else:
        # Should never reach here
        raise IOError('Need to set -r !!')

    # Run the setup
    ps.run(setup_only=True, sort_dir=sort_dir, write_bkg_pairs=args['background'])

    return ps.fitstbl, (ps, output_path)

def write_setup(args, context):
    ps, output_path = context
    # Use PypeItMetaData to write the complete PypeIt file
    pypeit_file = os.path.join(output_path, '{0}.pypeit'.format(args['spectrograph']))
    config_list = [item.strip() for item in args['cfg_split'].split(',')]
    
    return ps.fitstbl.write_pypeit(pypeit_file, cfg_lines=ps.user_cfg, write_bkg_pairs=args['background'],
                                configs=config_list)

def redux(args):
    splitnm = os.path.splitext(args['pypeit_file'])
    if splitnm[1] != '.pypeit':
        msgs.error("Bad extension for PypeIt reduction file."+msgs.newline()+".pypeit is required")
    logname = splitnm[0] + ".log"
    
    pypeIt = pypeit.PypeIt(args['pypeit_file'], verbosity=2,
                           reuse_masters=True, #~args['do_not_reuse_masters'],
                           overwrite=True, #args['overwrite'],
                           redux_path=None, #args['redux_path'], # rethink these
                           calib_only=False, #args['calib_only'],
                           logname=logname, show=False) #args['show'])
    pypeIt.reduce_all()
    msgs.info('Data reduction complete')

    msgs.info('Generating QA HTML')
    pypeIt.build_qa()