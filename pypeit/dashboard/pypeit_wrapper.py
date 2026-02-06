# I will make a wrapper class for pypeit and in this it will inherit all the functions in the pypeit class 
# but when I want to run it I will overide the run function and have it send me the state the program is at

from pypeit.pypeit import *
from PyQt6.QtCore import QObject, pyqtSignal # might need something like zmq to do this instead
from multiprocessing import Process

# instead of pyqt6 signalling I will have to use sockets to publish what step Pypeit is on

class PypeItWrapper(PypeIt):
    def __init__(self, pypeit_file, msg_queue, overwrite=True, reuse_calibs=False, show=False, redux_path=None, calib_only=False ):
        super().__init__(pypeit_file, overwrite, reuse_calibs, show, redux_path, calib_only)
        self.msg_queue = msg_queue
        self.calib_only = calib_only

        # -----------------------------------------------------------
        self.msg_queue.put(("Step","initializing pypeit"))
        # -----------------------------------------------------------

    def calib_all(self):

        # -----------------------------------------------------------
        self.msg_queue.put(("Step","Calibrations"))
        # -----------------------------------------------------------

        """
        Process all calibration frames.

        Provides an avenue to process the calibrations for a dataset 
        without (or omitting) any science/standard frames.

        I COPIED THESE STRAIGHT FROM PYPEIT SO IF ANYTHING CHANGES THERE IT WILL HAVE TO 
        MANUALLY BE CHANGED HERE
        """
        self.tstart = time.perf_counter()

        # Frame indices
        for calib_ID in self.fitstbl.calib_groups:

            # ------------------------------------------------------
            self.msg_queue.put(("Calibration ID",f'{calib_ID}'))
            # ------------------------------------------------------

            # Find all the frames in this calibration group
            in_grp = self.fitstbl.find_calib_group(calib_ID)
            if not any(in_grp):
                continue
            # Find the detectors to reduce
            detectors = self.spectrograph.select_detectors(subset=self.par['rdx']['detnum'] if self.par['rdx']['slitspatnum'] is None 
                                              else self.par['rdx']['slitspatnum'])
            log.info(f'Detectors to work on: {detectors}')

            # Loop on Detectors
            for self.det in detectors:

                # ------------------------------------------------------
                self.msg_queue.put(("Detector",f'{self.det}'))
                # ------------------------------------------------------

                log.info(f'Working on detector {self.det}')

                caliBrate = pypeit_steps.calib_one(self.spectrograph, self.fitstbl, self.par,
                                       self.det, calib_ID, self.calibrations_path)
                                       

        # Finish
        self.print_end_time()

    def reduce_all(self):
        
        # ------------------------------------------------------------
        self.msg_queue.put(("Step","Reductions")) 
        # -------------------------------------------------------------

        """
        Main driver of the end-to-end reduction

        Calibration and extraction via a series of calls to
        :func:`reduce_exposure`.

        """
        # Validate the parameter set
        self.par.validate_keys(required=['rdx', 'calibrations', 'scienceframe', 'reduce',
                                         'flexure'])
        self.tstart = time.perf_counter()

        # ############################################################################
        # Standard Star(s) Loop
        # ############################################################################
        # Iterate over each calibration group and reduce the standards
        for calib_ID in self.fitstbl.calib_groups:


            # ------------------------------------------------------------
            self.msg_queue.put(("Calibration ID",f"{calib_ID}")) 
            # -------------------------------------------------------------

            reduce_calibID(self.spectrograph, self.par, self.fitstbl,
                           calib_ID, self.calibrations_path,
                           reduce_standard=True, overwrite=self.overwrite,
                           show=self.show, 
                           run_state=self.run_state,
                           reuse_calibs=self.reuse_calibs)

        # ############################################################################
        # Science Frame(s) Loop
        # ############################################################################
        # Iterate over each calibration group again and reduce the science frames
        for calib_ID in self.fitstbl.calib_groups:
            reduce_calibID(self.spectrograph, self.par, self.fitstbl,
                           calib_ID, self.calibrations_path,
                           reduce_standard=False, overwrite=self.overwrite,
                           show=self.show, run_state=self.run_state,
                           reuse_calibs=self.reuse_calibs)
            log.info(f'Finished calibration group {calib_ID}')

        # Finish
        self.print_end_time()

    def run(self):
        # there is no logging for this but I have mostly copied the run function from run_pypeit.py
        # I will probably have to copy the logging but we shall see
        if self.calib_only:
            self.calib_all()
        else:
            self.reduce_all()

        self.build_qa()

        return 0



