from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import time
from abc import ABCMeta
import os
import datetime

from pypeit import msgs
from pypeit.core import fsort
from pypeit.core import qa
from pypeit import arms
from pypeit.par.util import make_pypeit_file, parse_pypeit_file

from pypeit import pypeitsetup
from pypeit.spectrographs import keck_lris
from pypeit.scripts import run_pypeit

from pypeit import debugger


class PypeIt(object):
    """
    This class is designed to run PypeIt

    .. todo::
        Improve docstring...

    Parameters
    ----------

    Attributes
    ----------
        pypeit_file (:obj:`str`):
            Name of the pypeit file to read.  Pypit files have a specific
            set of valid formats. A description can be found `here`_
            (include doc link).

    Inherited Attributes
    --------------------
    """
    __metaclass__ = ABCMeta

    def __init__(self, spectrograph, verbosity=2, overwrite=True, logname=None):

        # Init
        self.spectrograph = spectrograph
        self.verbosity = verbosity
        self.overwrite = overwrite

        # Internals
        self.pypeit_file = ''
        self.logname = logname


    def build_setup_files(self, files_root, redux_path=None):
        pargs, _ = self.make_setup_file(files_root, redux_path=redux_path)

        self._setup(setup_only=True, calibration_check=False)

    def build_custom_pypeitfile(self):

        msgs.reset(verbosity=2)



    def make_setup_file(self, files_root, extension='.fits', redux_path=None,
                        overwrite=False):

        # setup_files dir
        redux_path = os.getcwd() if redux_path is None else redux_path
        outdir = os.path.join(redux_path, 'setup_files')
        msgs.info('Setup files will be written to: {0}'.format(outdir))
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        # Generate a dummy .pypeit file
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        root = self.spectrograph.spectrograph+'_'+date
        pypeit_file = outdir+'/'+root+'.pypeit'

        # Generate
        dfname = "{:s}*{:s}*".format(files_root, extension)
        # configuration lines
        cfg_lines = ['[rdx]']
        cfg_lines += ['    spectrograph = {0}'.format(self.spectrograph.spectrograph)]
        cfg_lines += ['    sortroot = {0}'.format(root)]
        make_pypeit_file(pypeit_file, self.spectrograph, [dfname], cfg_lines=cfg_lines, setup_mode=True)
        msgs.info('Wrote template pypeit file: {0}'.format(pypeit_file))

        # Parser
        pinp = [pypeit_file, '-p', '-s {0}'.format(root) ]
        if overwrite:
            pinp += ['-o']
        pargs = run_pypeit.parser(pinp)
        self.sorted_file = pypeit_file.replace('.pypeit', '.sorted')
        # Return
        return pargs, self.sorted_file

    def msgs_reset(self):
        # Reset the global logger
        msgs.reset(log=self.logname, verbosity=self.verbosity)
        msgs.pypeit_file = self.pypeit_file

    def run(self, quick=False):
        """
        Execute PYPIT.

        .. todo::
            - More description in docstring
            - Allow the user to just provide a list of files or always
              require a pypeit file?

        Args:
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
            raise NotImplementedError('Quick version of pypeit is not yet implemented.')

        # Msgs
        self.msgs_reset()

        # Setup
        self.setup()

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
        fsort.make_dirs(spectrograph.spectrograph, par['calibrations']['caldir'],
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
            qa.gen_mf_html(pypeit_file)
            qa.gen_exp_html()
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

    def _setup(self, setup_only=False, calibration_check=False, use_header_frametype=False, sort_dir=None):
        """

        Args:
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

        Returns:

        """

        # Msgs
        self.msgs_reset()

        # Record the starting time
        tstart = time.time()

        # Perform the setup
        self.setup = pypeitsetup.PypeItSetup.from_pypeit_file(self.pypeit_file)
        par, _, self.fitstbl, self.setup_dict = self.setup.run(setup_only=setup_only,
                                                           calibration_check=calibration_check,
                                                           use_header_frametype=use_header_frametype,
                                                           sort_dir=sort_dir)
        # Write the fits table
        self.setup.write_fitstbl()

class LRISb(PypeIt):
    def __init__(self):
        spectrograph = keck_lris.KeckLRISBSpectrograph()
        PypeIt.__init__(self, spectrograph)



