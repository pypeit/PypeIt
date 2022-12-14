import argparse

from pypeit.scripts import scriptbase
from pypeit.setup_gui.view import PypeItSetupView
from pypeit.setup_gui.model import PypeItSetupProxy
from pypeit.setup_gui.controller import SetupGUIController

        

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
        
        controller = SetupGUIController(args, qt_args)
        controller.start()

if __name__ == '__main__':
    SetupGUI.entry_point()
