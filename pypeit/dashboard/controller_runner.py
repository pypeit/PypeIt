import sys
from PyQt6.QtWidgets import QApplication
from pypeit.setup_gui.controller import SetupGUIController

app = QApplication(sys.argv)
verbosity = 1 # in the future I will store the verbosity in its own file and then read the verbosity from that file maybe
controller = SetupGUIController(app,verbosity)

controller.start()
