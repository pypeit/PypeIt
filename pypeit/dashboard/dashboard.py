import sys
import io
import re
import subprocess
from pathlib import Path
from qtpy import QtWidgets
import logging
from logging.handlers import QueueListener
from multiprocessing import Process, Queue
from contextlib import redirect_stdout

from PyQt6.QtCore import pyqtSignal
import qtpy
from qtpy.QtCore import QTimer, QSize, Qt, QMargins, QObject
from qtpy.QtGui import QIcon, QColor, QColorConstants, QPainter
from qtpy.QtWidgets import QApplication, QWidget, QHBoxLayout, QVBoxLayout, QPushButton, QGridLayout, QLabel, QProgressBar,QTabWidget, QListWidget, QAbstractItemView, QFileDialog, QTextEdit

# from pypeit.setup_gui.controller import start_gui
# from pypeit.scripts import setup
import pypeit
from pypeit.dashboard.pypeit_worker import PypeItWorker, check_pypeit_status 
# from pypeit.dashboard.step_tracker import message_listener

"""
TODO: give a meta view and specific (show meta step, what step of that step are we one, what step of that step are we on)
make another tab next to qa and science that shows the logs of run_pypeit

NOTE: args is argparse so the arguments you put in the terminal
"""

class FilledBackgroundWidget(QWidget):
    def __init__(self, color=None):
        super().__init__()
        self.color = color

    def paintEvent(self, event):
        super().paintEvent(event)
        if self.color is not None:
            painter = QPainter(self)
            painter.fillRect(event.rect(), QColor(self.color))
    

class ButtonWidget(FilledBackgroundWidget):
    """this widget is the top left corner"""
    def __init__(self):
        super().__init__()#color=QColorConstants.DarkBlue)

        #----------------- defining the widgets ----------------------
        self.open_setup_button = QPushButton()
        self.edit_setup_button = QPushButton()
        self.run_all_button = QPushButton()
        self.run_next_button = QPushButton()
        self.help_button = QPushButton()
        self.check_status_button = QPushButton()

        self.test_icons = [(QIcon.ThemeIcon.DocumentOpen, "Open Setup", self.open_setup_button),
                           (QIcon.ThemeIcon.InputKeyboard, "Import Setup", self.edit_setup_button),
                           (QIcon.ThemeIcon.MediaSeekForward, "Run All", self.run_all_button),
                           (QIcon.ThemeIcon.MediaSkipForward, "Run Next", self.run_next_button),
                           (QIcon.ThemeIcon.HelpFaq, "Help", self.help_button),
                           (QIcon.ThemeIcon.Computer, "Check Status", self.check_status_button),
                           ]

        layout=QVBoxLayout()

        for icon, text, widget in self.test_icons:
            widget.setStyleSheet(f"text-align:left;")

            widget.setText(text)
            widget.setIcon(QIcon.fromTheme(icon))
            widget.setIconSize(QSize(32,32))
            layout.addWidget(widget)
            

        self.setLayout(layout)
        self.layout().setContentsMargins(0, 0, 0, 0)
        self.setMaximumWidth(self.fontMetrics().averageCharWidth()*20)
        self.setMaximumHeight(self.fontMetrics().lineSpacing()*15)




class StatusWidget(FilledBackgroundWidget):
    """this widget is a collection of widgets at the top middle (from setup file to status bar)"""
    """value_style_sheet is something that makes a little box underneath it"""

    def __init__(self):
        super().__init__()#color=QColorConstants.DarkGreen)
        fm = self.fontMetrics()       
        h = fm.lineSpacing() * 8
        #self.setMinimumHeight(h)
        self.setMaximumHeight(h)
        w = fm.averageCharWidth() * 80
        self.setMaximumWidth(w)

        value_cm = QMargins(fm.averageCharWidth(),0,fm.averageCharWidth(),0)
        value_style_sheet = "background-color:rgb(80,80,80);"
        layout = QGridLayout()
        #l = QLabel(text="Status")
        #l.setStyleSheet("font: normal 36pt")
        #layout.addWidget(l, 0, 0, 1, 3,alignment=Qt.AlignmentFlag.AlignLeft)
        

        #---------------------- setup file group -------------------
        setup_file_label = QLabel(text="Setup File")
        layout.addWidget(setup_file_label,0,0,1,1)#,alignment=Qt.AlignmentFlag.AlignLeft)

        self.setup_file = QLabel(text="None")
        self.setup_file.setContentsMargins(value_cm)
        self.setup_file.setStyleSheet(value_style_sheet)
        layout.addWidget(self.setup_file,0,1,1,1)
        
        # --------------------- Calibration id group --------------
        calibration_id_label = QLabel(text="Calibration ID")
        layout.addWidget(calibration_id_label,1,0,1,1)

        self.calibration_id = QLabel(text="None")
        self.calibration_id.setContentsMargins(value_cm)
        self.calibration_id.setStyleSheet(value_style_sheet)
        layout.addWidget(self.calibration_id,1,1,1,1)#,alignment=Qt.AlignmentFlag.AlignLeft)
        
        # --------------------- Detector group ---------------------
        detector_label = QLabel(text="Detector")
        layout.addWidget(detector_label,2,0,1,1)

        self.detector = QLabel(text="None")
        self.detector.setContentsMargins(value_cm)
        self.detector.setStyleSheet(value_style_sheet)
        layout.addWidget(self.detector,2,1,1,1)
        
        # ------------------- Science file group ------------------
        science_file_label = QLabel(text="Science File")
        layout.addWidget(science_file_label,0,2,1,1)

        self.science_file = QLabel(text="None")
        self.science_file.setContentsMargins(value_cm)
        self.science_file.setStyleSheet(value_style_sheet)
        layout.addWidget(self.science_file,0,3,1,1)

        # ------------------- step group ----------------------
        step_label = QLabel(text="Step")
        layout.addWidget(step_label,1,2,1,1)

        self.meta_step = QLabel(text="None")
        self.meta_step.setContentsMargins(value_cm)
        self.meta_step.setStyleSheet(value_style_sheet)
        layout.addWidget(self.meta_step,1,3,1,1)

        # ------------------------ calibration step group ----------------
        calibration_step_label = QLabel(text="Calibration Step")
        layout.addWidget(calibration_step_label,2,2,1,1)

        self.calibration_step = QLabel(text="None")
        self.calibration_step.setContentsMargins(value_cm)
        self.calibration_step.setStyleSheet(value_style_sheet)
        layout.addWidget(self.calibration_step,2,3,1,1)

        # ------------------ progress bar ---------------------
        progress_bar = QProgressBar()
        progress_bar.setMaximum(100)
        progress_bar.setValue(33)
        progress_bar.setTextVisible(True)
        layout.addWidget(progress_bar,3,0,1,4)

        #layout.addWidget(SpacerWidget(rows=5,cols=40),3,4,3,4)
        layout.setVerticalSpacing(self.fontMetrics().lineSpacing())
        layout.setHorizontalSpacing(self.fontMetrics().averageCharWidth())
        self.setLayout(layout)
        cm = self.layout().contentsMargins()
        cm.setTop(0)
        #self.layout().setContentsMargins(cm)
        
    def update_setup_file(self,setup_file_path):
        # sets the setup_file_label to have the setup file path next to it that is updated
        self.setup_file.setText(str(setup_file_path))
    def update_calibration_id(self,step):
        self.calibration_id.setText(step)
    def update_detector_step(self,step):
        self.detector.setText(step)
    def update_science_file(self,science_file):
        self.science_file.setText(str(science_file))
    def update_meta_step(self,step:str):
        self.meta_step.setText(step)
    def update_calibration_step(self,step):
        self.calibration_step.setText(step)
    def update_progress_bar(self,update):
        pass # I don't really know how much I will need this but we shall see


        

# special QlistWidget for easy adding of items
class logs_view_widget(QTextEdit):
    def __init__(self):
        super().__init__()
        self.setReadOnly(True)
        self.setLineWrapMode(QtWidgets.QTextEdit.LineWrapMode.NoWrap)

    def update_logs(self,message):
        self.append(message)


class DashboardWidget(FilledBackgroundWidget):
    """this widget is the bottom widget"""
    def __init__(self):
        super().__init__()

        layout = QVBoxLayout()

        # ------------- definitions --------------
        self.status_widget = StatusWidget()
        self.calibration_widget = QTextEdit() # should eventually change this to a table
        self.qa_widget = QListWidget()
        self.science_widget = QListWidget()
        self.logs_widget = logs_view_widget()

        # -------------- main layout -------------
        layout.addWidget(self.status_widget)
        tab_widget = QTabWidget()
        tab_widget.addTab(self.qa_widget,"QA") # Quality analysis
        tab_widget.addTab(self.calibration_widget,"Calibrations")
        tab_widget.addTab(self.science_widget,"Science")
        tab_widget.addTab(self.logs_widget,"Logs")

        layout.addWidget(tab_widget, 3)
        self.setLayout(layout)


def parse_pypeit_setup_file(file_path):
    """
    takes in a file path for a pypeit_setup file

    args: file path of the pypeit_setup file

    returns: tuple of (spectrograph, raw_path)
    """
    with open(file_path,"r") as setup:
        # contents = setup.readlines()
        spectrograph = None
        raw_path = None
        science_file = None
        spectrograph_search_string = r'^\s*spectrograph\s*=\s*(\s+)' # this should give the spectrograph after the word spectrograph = 
        file_path_search_string = r'^\s*#?\s*path\s+(.+)' # thi should give the path after the word path
        other_file_pattern = re.compile(r'^\|\s*(\S+\.fits)\s*\|\s*([^|]+?)\s*(?:\||$)')


        while spectrograph == None:
            line = setup.readline() # I am pretty sure this will go through each line instead of the same line
            spectrograph = re.search(spectrograph_search_string,line)

        while raw_path == None:
            line = setup.readline()
            raw_path = re.search(file_path_search_string,line)

        files = [other_file_pattern.match(line) for line in setup.readlines()]
        files = [(x.group(1),x.group(2).strip()) for x in files if x]

        # get the science file, this could cause an error so put in a try except 
        try:
            science_file = [file[0] for file in files if file[1] == "science"][0]
        except:
            pass

        spectrograph = spectrograph.group(1)
        raw_path = raw_path.group(1)
        # now its time to get the files from the pypeit setup file

    return spectrograph,raw_path,files,science_file

# this is to handle the logging from pypeit. will be redirected to the dashboard for it to display

class MainWindow(QWidget):
    
    def __init__(self):
        super().__init__()

        layout = QHBoxLayout()
        self.setup_widget = ButtonWidget()
        self.dashboard_widget = DashboardWidget()
        layout.addWidget(self.setup_widget,alignment=Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self.dashboard_widget,stretch=3)

        self.setup_file_path = None

        # -------- connections ---------
        self.setup_widget.open_setup_button.clicked.connect(self.start_controller)
        self.setup_widget.edit_setup_button.clicked.connect(self.import_setup_file)
        self.setup_widget.check_status_button.clicked.connect(self.check_status)
        # run_all_button clicked is handled in main. soon I will maybe make all of 
        # the connections to be handled in main but who knows


        self.setLayout(layout)

    def start_controller(self):
        subprocess.Popen(["pypeit_setup","--gui"]) # starting the controller runnner file

    def check_status(self):
        # need to have a seperate thing for when pypeit is running or not
        if self.setup_file_path != None:
            f = io.StringIO()
            with redirect_stdout(f):
                check_pypeit_status(self.setup_file_path)
            captured_output = f.getvalue()
            self.dashboard_widget.calibration_widget.append(captured_output)



 
    def import_setup_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select a PypeIt setup file",
            "",
            "PypeIt (*.pypeit)" 
        )

        if file_path: 
            # maybe add a fuction to this class that just updates everything by calling other functions for the class
            # I should update this so that it is something easier understood and doesn't return too many things at once 
            spectrograph,raw_path,calibration_files,science_file = parse_pypeit_setup_file(file_path) 

            file_name = Path(file_path).name
            self.dashboard_widget.status_widget.update_setup_file(file_name)
            self.update_science_file(science_file)
            self.setup_file_path = file_path

    def update_science_file(self,science_file):
        self.dashboard_widget.status_widget.update_science_file(science_file) 

    def update_logs(self, line):
        self.dashboard_widget.logs_widget.update_logs(line)

    def update_qa_tab(self,directory_path):
        items = [str(item) for item in directory_path.iterdir()]
        self.dashboard_widget.qa_widget.addItems(items)
    def update_science_tab(self,directory_path):
        items = [str(item) for item in directory_path.iterdir()]
        self.dashboard_widget.science_widget.addItems(items)


# ------------------------------------------------------------------------
# these are classes and functions for multiprocessing and threading
class QtLogHandler(QObject, logging.Handler):

    log_signal = pyqtSignal(str)
    step_signal = pyqtSignal(str)

    def __init__(self):
        QObject.__init__(self)
        logging.Handler.__init__(self)

    def emit(self, record):
        msg = self.format(record)
        
        if record.levelname == "STEP":
            self.step_signal.emit(msg)
        else:
            self.log_signal.emit(msg)

def start_pypeit_process(file_path,log_queue):
    p = Process(target=PypeItWorker,args=(f'{file_path}',log_queue),daemon=True)
    p.start()
    return p

def step_listener(step):
    # print(step)
    pass

# --------------------------------------------
def main():
        # Note QT expects the program name as arg 0
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    app = QApplication(sys.argv)

    # Setup application/window icon
    iconPath = Path(__file__).parent.parent / "setup_gui/images/window_icon.png"
    if not iconPath.exists():
        num = 1
    else:
        app.setWindowIcon(QIcon(str(iconPath)))
    

    defaultFont = app.font()
    if defaultFont.pointSizeF() < 18.0:
        defaultFont.setPointSize(18)
        app.setFont(defaultFont)

    main_window = MainWindow()
    main_window.setWindowTitle(main_window.tr("PypeIt Dashboard"))
    main_window.resize(1650,900)
    main_window.show()

    # ---------------- LOGGING ----------------
    log_queue = Queue()

    qt_handler = QtLogHandler()
    qt_handler.setFormatter(logging.Formatter(
        " %(levelname)s | %(message)s"
    ))

    qt_handler.log_signal.connect(main_window.update_logs)

    log_listener = QueueListener(log_queue, qt_handler)
    log_listener.start()
    # ----------------- STEP TRACKER ----------


    main_window.setup_widget.run_all_button.clicked.connect(
            lambda: start_pypeit_process(
                main_window.setup_file_path,
                log_queue
            )
        )
    # --------------------- this is for the SetupGUIController ----------------         

    # QT runs it's event loop in C, so the python signal handling mechanism
    # is never called, or it's only called after you give focus to the
    # window. To make Ctrl+C handling work immediately in a way that still 
    # calls the PypeIt CTRL+C handler, we set a timer to run every 500ms in the
    # python interpreter, which will allow the python signal handling
    # code to it.
        
    # This trck was brought to you by this stack exchange thread:
    # https://stackoverflow.com/questions/4938723/what-is-the-correct-way-to-make-my-pyqt-application-quit-when-killed-from-the-co
    timer = QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)

    # Start the applications event loop
    app.exec()


if __name__ == '__main__':
    sys.exit(main())
