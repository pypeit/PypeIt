# I will make a wrapper class for pypeit and in this it will inherit all the functions in the pypeit class 
# but when I want to run it I will overide the run function and have it send me the state the program is at

from pypeit.pypeit import PypeIt
from PyQt6.QtCore import QObject, pyqtSignal # might need something like zmq to do this instead
from multiprocessing import Process

# instead of pyqt6 signalling I will have to use sockets to publish what step Pypeit is on
class SignalSender(QObject):
    calib_started = pyqtSignal(str)
    reduce_started = pyqtSignal(str)

class PypeItWrapper(PypeIt):
    def __init__(self, pypeit_file, verbosity=2, overwrite=True, reuse_calibs=False, logname=None, show=False, redux_path=None, calib_only=False):
        super().__init__(pypeit_file, verbosity, overwrite, reuse_calibs, logname, show, redux_path, calib_only)
    # simple class that adds a signal to when calib_all and reduce_all is called. 
    # later hopefully will be able to signal the individual steps
    def calib_all(self):
        # do something that will say that we started calibrations
        return super().calib_all()

    def reduce_all(self):
        # say something that will say we started reduction
        return super().reduce_all()

    def run(self):
        # there is no logging for this but I have mostly copied the run function from run_pypeit.py
        # I will probably have to copy the logging but we shall see
        if self.calib_only:
            self.calib_all()
        else:
            self.reduce_all()

        self.build_qa()

        return 0



