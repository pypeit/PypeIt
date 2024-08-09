"""
An interactive GUI for creaing pypeit input files.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import sys
import os
from qtpy.QtWidgets import QApplication

from pypeit.scripts import scriptbase
from pypeit.setup_gui.controller import SetupGUIController
from pypeit.spectrographs import available_spectrographs
class SetupGUI(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description="Interactive GUI for creating and editing PypeIt input files. "
                                                "Additional Qt arguments can also be used. See https://doc.qt.io/qt-5/qapplication.html#QApplication", 
                                    width=width)
        parser.add_argument('-s', '--spectrograph', default=None, type=str,
                            help='A valid spectrograph identifier: {0}'.format(
                                    ', '.join(available_spectrographs)))
        parser.add_argument('-r', '--root', default=[], type=str,nargs='+',
                            help='Root to search for data files.  You can provide the top-level '
                                 'directory  (e.g., /data/Kast) or the search string up through '
                                 'the wildcard (.e.g, /data/Kast/b).  Use the --extension option '
                                 'to set the types of files to search for.  Default is the '
                                 'current working directory.')
        parser.add_argument('-e', '--extension', default='.fits',
                            help='File extension; compression indicators (e.g. .gz) not required.')
        parser.add_argument('-l', '--logfile', type=str, default=None, 
                            help="Write the PypeIt logs to the given file. If the file exists it will be renamed.")
        parser.add_argument('-v', '--verbosity', type=int, default=2,
                            help='Verbosity level between 0 [none] and 2 [all]. Default: 2.')
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
        # Set the Qt Arguments. Note QT expects the program name as arg 0
        qt_args = [sys.argv[0]] + combined_args[1]
        app = QApplication(qt_args)
        controller = SetupGUIController(args)
        controller.start(app)
