import signal
import sys
import argparse
import datetime

from pypeit.scripts import scriptbase
from pypeit.setup_gui.view import PypeItSetupView
from pypeit.setup_gui.model import PypeItSetupModel
from pypeit.setup_gui.controller import PypeItSetupController

from qtpy.QtCore import QCoreApplication
from qtpy.QtWidgets import QApplication
        

class SetupGUI(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = argparse.ArgumentParser("Interactive GUI for creating and editing PypeIt input files.",
                                         epilog="Additional Qt arguments can also be used. See https://doc.qt.io/qt-5/qapplication.html#QApplication")
        parser.add_argument('-v', '--verbosity', type=int, default=2,
                            help='Verbosity level between 0 [none] and 2 [all]. Default: 2. '
                                 'Level 2 writes a log with filename setup_gui_YYYYMMDD-HHMM.log')
        return parser

    @classmethod
    def parse_args(cls, options=None):
        """
        Parse the command-line arguments.
        """
        parser = cls.get_parser()
        # Use parse_known args so we can pass the remainder over to Qt
        return  parser.parse_known_args() if options is None else parser.parse_known_args(options)        

    @staticmethod
    def main(combined_args):

        args = combined_args[0]
        qt_args = combined_args[1]

        app = QApplication(qt_args)

        QCoreApplication.setOrganizationName("PypeIt")
        QCoreApplication.setApplicationName("SetupGUI")
        QCoreApplication.setOrganizationDomain("pypeit.readthedocs.io")
        
        timestamp = datetime.datetime.utcnow().strftime("%Y%m%d-%H%M")
        logname = f"setup_gui_{timestamp}.log"
        model = PypeItSetupModel(log_file=logname, verbosity=args.verbosity)

        main_window = PypeItSetupView()
        main_window.resize(1650,900)

        controller = PypeItSetupController(model, main_window)
        main_window.show()

        # QT Gobbles up the Python Ctrl+C handling, so the default PypeIt Ctrl+C handler won't
        # be called. So we reset it to the OS default
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        sys.exit(app.exec_())

if __name__ == '__main__':
    SetupGUI.entry_point()
