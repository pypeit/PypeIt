
# this is the worker process and it will process the logs



import logging
from logging.handlers import QueueHandler
from pypeit import log
from pypeit.scripts.run_pypeit import RunPypeIt
from pypeit import pypeit

def run_pypeit_main_wrapper(args):
    pypeIt = pypeit.PypeIt(
        args.pypeit_file, reuse_calibs=args.reuse_calibs, overwrite=args.overwrite,
        redux_path=args.redux_path, calib_only=args.calib_only, show=args.show
    )

    if args.calib_only:
        pypeIt.calib_all()
    else:
        pypeIt.reduce_all()
    log.info('Data reduction complete')

    # QA HTML
    log.info('Generating QA HTML')
    pypeIt.build_qa()

    return 0

def PypeItWorker(setup_file_path,log_queue):


    parser = RunPypeIt.get_parser()

    args = parser.parse_args([
        setup_file_path,
    ])

    # Run the worker
    # RunPypeIt.main(args)
    return run_pypeit_main_wrapper(args)

def check_pypeit_status(setup_file_path):

    parser = RunPypeIt.get_parser()

    args = parser.parse_args([
        setup_file_path,
    ])

    pypeIt = pypeit.PypeIt(
            args.pypeit_file,
            reuse_calibs=True,
            calib_only=True)

    pypeIt.calib_all(status_only=True, reload_only=True)

    # Write state to JSON
    pypeIt.run_state.write()

    # Pretty-print the state
    pypeIt.run_state.print_status()
