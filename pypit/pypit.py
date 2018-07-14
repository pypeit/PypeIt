from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import time

from pypit import msgs
# TODO: (KBW) Can archeck code be put in pypit/__init__.py ?
from pypit import archeck  # THIS IMPORT DOES THE CHECKING.  KEEP IT
from pypit.core import arsort
from pypit import arqa
from pypit import arms

from pypit import pypitsetup

from pypit import ardebug

def PYPIT(pypit_file, setup_only=False, calibration_check=False, use_header_frametype=False,
          sort_dir=None, debug=None, quick=False, ncpus=1, overwrite=False, verbosity=1,
          use_masters=False, logname=None):
    """
    Execute PYPIT.

    .. todo::
        - More description in docstring
        - Allow the user to just provide a list of files or always
          require a pypit file?

    Args:
        pypit_file (:obj:`str`):
            Name of the pypit file to read.  Pypit files have a specific
            set of valid formats. A description can be found `here`_
            (include doc link).
        setup_only (bool):
            Only this setup will be performed.  Pypit is expected to
            execute in a way that ends after this class is fully
            instantiated such that the user can inspect the results
            before proceeding.  This has the effect of providing more
            output describing the success of the setup and how to
            proceed, and provides warnings (instead of errors) for
            issues that may cause the reduction itself to fail.
        calibration_check (bool):
            Only check that the calibration frames are appropriately
            setup and exist on disk.  Pypit is expected to execute in a
            way that ends after this class is fully instantiated such
            that the user can inspect the results before proceeding. 
        use_header_frametype (bool):
            Allow setup to use the frame types drawn from the file
            headers using the instrument specific keywords.
        sort_dir (str):
            The directory to put the '.sorted' file.
        debug (:obj:`dict`, optional):
            Debugging dictionary.  TODO: More description.
        quick (:obj:`bool`, optional):
            Perform a quick version of the reduction.  NOT IMPLEMENTED.
        ncpus (:obj:`int`, optional):
            The number of cpus to use.  NOT IMPLEMENTED.
        overwrite (:obj:`bool`, optional):
            Flag to overwrite any existing files/directories.
        verbosity (:obj:`int`, optional):
            Verbosity level of system output.  Can be::
                - 0: No output
                - 1: Minimal output (default)
                - 2: All output
        use_masters (:obj:`bool`, optional):
            Use the master frames if available (same as setting
            par['calibrations']['masters'] = 'reuse'.  NOT IMPLEMENTED.
        logname (:obj:`str`, optional):
          The name of an ascii log file with the details of the
          reduction

    Returns:
        int: The status of the reduction::
            - 0: Reductions successful
            - 1: Setup successful (when `setup_only=True`)
            - 2: Calibration check successful (when `calibration_check=True`)
    """
    if quick:
        raise NotImplementedError('Quick version of pypit is not yet implemented.')

    # Init logger
    if debug is None:
        debug = ardebug.init()

    # Reset the global logger
    msgs.reset(log=logname, debug=debug, verbosity=verbosity)
    msgs.pypit_file = pypit_file

    # Record the starting time
    tstart = time.time()

    # Perform the setup
    setup = pypitsetup.PypitSetup.from_pypit_file(pypit_file)
    par, spectrograph, fitstbl, setup_dict = setup.run(setup_only=setup_only,
                                                       calibration_check=calibration_check,
                                                       use_header_frametype=use_header_frametype,
                                                       sort_dir=sort_dir)
    # Write the fits table
    setup.write_fitstbl()

    # Exit if finished
    if setup_only:
        msgs.info('Setup complete')
        return 1
    if calibration_check:
        msgs.info('Calcheck complete')
        return 2

    # Make the output directories
    # TODO: Do we want the code to interactively ask for a new
    # directory?  I think it would be better if it just faulted when a
    # directory/file exists and overwrite is False.
    arsort.make_dirs(spectrograph.spectrograph, par['calibrations']['caldir'],
                     par['rdx']['scidir'], par['rdx']['qadir'], overwrite=overwrite)

    # Just do it (sponsored by Nike)
    if par['rdx']['pipeline'] == 'ARMS':
        msgs.info('Data reduction will be performed using PYPIT-ARMS')
        #status = arms.ARMS(fitstbl, setup_dict, sciexp=sciexp)
        status = arms.ARMS(fitstbl, setup_dict, par=par, spectrograph=spectrograph)
    elif par['rdx']['pipeline'] == 'ARMED':
        import pdb; pdb.set_trace()
        msgs.error('ARMED is currently broken.')
        msgs.info('Data reduction will be performed using PYPIT-ARMED')
        status = armed.ARMED(fitstbl)
    else:
        msgs.error('Unrecognized pipeline!')

    # Check for successful reduction
    if status == 0:
        msgs.info('Data reduction complete')
        # QA HTML
        msgs.info('Generating QA HTML')
        arqa.gen_mf_html(pypit_file)
        arqa.gen_exp_html()
    else:
        msgs.error('Data reduction failed with status ID {0:d}'.format(status))

    # Capture the end time and print it to user
    tend = time.time()
    codetime = tend-tstart
    if codetime < 60.0:
        msgs.info('Data reduction execution time: {0:.2f}s'.format(codetime))
    elif codetime/60.0 < 60.0:
        mns = int(codetime/60.0)
        scs = codetime - 60.0*mns
        msgs.info('Data reduction execution time: {0:d}m {1:.2f}s'.format(mns, scs))
    else:
        hrs = int(codetime/3600.0)
        mns = int(60.0*(codetime/3600.0 - hrs))
        scs = codetime - 60.0*mns - 3600.0*hrs
        msgs.info('Data reduction execution time: {0:d}h {1:d}m {2:.2f}s'.format(hrs, mns, scs))

    return status

