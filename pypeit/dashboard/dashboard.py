import sys
import subprocess
import re
from pathlib import Path
from pypeit import msgs
import zmq

from qtpy.QtCore import QTimer, QSize, Qt, QMargins
from qtpy.QtGui import QIcon, QColor, QColorConstants, QPainter
from qtpy.QtWidgets import QApplication, QWidget, QHBoxLayout, QVBoxLayout, QPushButton, QGridLayout, QLabel, QProgressBar,QTabWidget, QListWidget, QAbstractItemView, QFileDialog
import qtpy

# from pypeit.setup_gui.controller import start_gui
# from pypeit.scripts import setup
import pypeit
from pypeit.dashboard.capture_logs import PypeitWorker

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

        self.test_icons = [(QIcon.ThemeIcon.DocumentOpen, "Open Setup", self.open_setup_button),
                           (QIcon.ThemeIcon.InputKeyboard, "Import Setup", self.edit_setup_button),
                           (QIcon.ThemeIcon.MediaSeekForward, "Run All", self.run_all_button),
                           (QIcon.ThemeIcon.MediaSkipForward, "Run Next", self.run_next_button),
                           (QIcon.ThemeIcon.HelpFaq, "Help", self.help_button),
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
        msgs.info(f"status height {h} width{w}")

        value_cm = QMargins(fm.averageCharWidth(),0,fm.averageCharWidth(),0)
        value_style_sheet = "background-color:rgb(80,80,80);"
        layout = QGridLayout()
        #l = QLabel(text="Status")
        #l.setStyleSheet("font: normal 36pt")
        #layout.addWidget(l, 0, 0, 1, 3,alignment=Qt.AlignmentFlag.AlignLeft)
        

        #---------------------- setup file group -------------------
        setup_file_label = QLabel(text="Setup File")
        layout.addWidget(setup_file_label,0,0,1,1)#,alignment=Qt.AlignmentFlag.AlignLeft)

        self.setup_file = QLabel(text="keck_deimos_830g_m_8500.pypeit")
        self.setup_file.setContentsMargins(value_cm)
        self.setup_file.setStyleSheet(value_style_sheet)
        layout.addWidget(self.setup_file,0,1,1,1)
        
        # --------------------- Calibration id group --------------
        calibration_id_label = QLabel(text="Calibration ID")
        layout.addWidget(calibration_id_label,1,0,1,1)

        self.calibration_id = QLabel(text="0")
        self.calibration_id.setContentsMargins(value_cm)
        self.calibration_id.setStyleSheet(value_style_sheet)
        layout.addWidget(self.calibration_id,1,1,1,1)#,alignment=Qt.AlignmentFlag.AlignLeft)
        
        # --------------------- Detector group ---------------------
        detector_label = QLabel(text="Detector")
        layout.addWidget(detector_label,2,0,1,1)

        self.detector = QLabel(text="3")
        self.detector.setContentsMargins(value_cm)
        self.detector.setStyleSheet(value_style_sheet)
        layout.addWidget(self.detector,2,1,1,1)
        
        # ------------------- Science file group ------------------
        science_file_label = QLabel(text="Science File")
        layout.addWidget(science_file_label,0,2,1,1)

        self.science_file = QLabel(text="DE.20100913.22358.fits.gz")
        self.science_file.setContentsMargins(value_cm)
        self.science_file.setStyleSheet(value_style_sheet)
        layout.addWidget(self.science_file,0,3,1,1)

        # ------------------- step group ----------------------
        step_label = QLabel(text="Step")
        layout.addWidget(step_label,1,2,1,1)

        self.meta_step = QLabel(text="Calibrations")
        self.meta_step.setContentsMargins(value_cm)
        self.meta_step.setStyleSheet(value_style_sheet)
        layout.addWidget(self.meta_step,1,3,1,1)

        # ------------------------ calibration step group ----------------
        calibration_step_label = QLabel(text="Calibration Step")
        layout.addWidget(calibration_step_label,2,2,1,1)

        self.calibration_step = QLabel(text="Tilts")
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
    def update_calibration_id(self,input):
        pass
    def update_detector_step(self,step):
        pass
    def update_science_file(self,science_file):
        self.science_file.setText(str(science_file))
    def update_meta_step(self,step):
        pass
    def update_calibration_step(self,step):
        pass
    def update_progress_bar(self,update):
        pass # I don't really know how much I will need this but we shall see


        
class FileListWidget(QListWidget):
    def __init__(self):
        super().__init__()
        self.addItems(["Arc_A_0_DET03.fits",
                       "Arc_A_0_DET07.fits",
                       "Edges_A_0_DET03.fits.gz",
                       "Edges_A_0_DET07.fits.gz",
                       "Slits_A_0_DET03.fits.gz",
                       "Slits_A_0_DET07.fits.gz",
                       "Tiltimg_A_0_DET03.fits",
                       "Tiltimg_A_0_DET07.fits",
                       "Tilts_A_0_DET03.fits",
                       "WaveCalib_A_0_DET03.fits,"
                       "WaveCalib_A_0_DET07.fits",
                       ])
        self.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)


class DashboardWidget(FilledBackgroundWidget):
    """this widget is the bottom widget"""
    def __init__(self):
        super().__init__()

        layout = QVBoxLayout()
        self.status_widget = StatusWidget()
        self.file_list_widget = FileListWidget()
        self.logs_widget = QListWidget()
        layout.addWidget(self.status_widget)
        tab_widget = QTabWidget()
        tab_widget.addTab(FilledBackgroundWidget(color=QColorConstants.Red),"QA") # Quality analysis
        tab_widget.addTab(self.file_list_widget,"Calibrations")
        tab_widget.addTab(FilledBackgroundWidget(color=QColorConstants.DarkGreen),"Science")
        tab_widget.addTab(self.logs_widget,"Logs")

        layout.addWidget(tab_widget, 3)
        self.setLayout(layout)



# what does each button do. 
# the open setup button opens the setup at the current stage/ state
# edit setup allows you to enter a different setup file
# run all starts running the tasks and run does this step by step
# help does something

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
        self.setup_widget.run_all_button.clicked.connect(self.run_all)


        self.setLayout(layout)

    def start_controller(self):
        subprocess.Popen(["pypeit_setup","--gui"]) # starting the controller runnner file

    def run_all(self):
        """
        how to get what the current step is that is happening in pypeit. 
        there is a pypeit_steps thing that runs_steps apparently and it also outputs what step it is on in logging
        I will regex the logs and when it finds pypeit_steps() or something similar, it will update the current step
        to be whatever comes after. 
        """

        # will need better checking in the future but for now will say if file path is not none
        command = ["run_pypeit",f"{self.setup_file_path}"]
        if self.setup_file_path != None:
            subprocess.Popen(command)
            worker = PypeitWorker(self.setup_file_path)
            worker.line_received.connect(self.update_logs)

            worker.run()
        else:
            # I will do something here that is like you don't have a setup file imported
            pass
 
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
            self.update_file_list_widget(calibration_files)
            self.update_science_file(science_file)
            self.setup_file_path = file_path


    def update_file_list_widget(self,file_list):
        # This currently will be input as the first element in each tuple will be the file
        new_list = [x[0] for x in file_list]
        self.dashboard_widget.file_list_widget.clear()
        self.dashboard_widget.file_list_widget.addItems(new_list)

    def update_science_file(self,science_file):
        self.dashboard_widget.status_widget.update_science_file(science_file) 

    def update_logs(self, line):
        self.dashboard_widget.logs_widget.addItems([line])


def main():
        # Note QT expects the program name as arg 0
    app = QApplication(sys.argv)

    # Setup application/window icon
    iconPath = Path(__file__).parent.parent / "setup_gui/images/window_icon.png"
    if not iconPath.exists():
        msgs.info("Icon path does not exist")
    else:
        app.setWindowIcon(QIcon(str(iconPath)))
    
    msgs.reset(verbosity=2, log="dashboard.log", log_to_stderr=True)
    msgs.info(f"QT Version: {qtpy.QT_VERSION}")
    msgs.info(f"PySide version: {qtpy.PYSIDE_VERSION}")
    msgs.info(f"PyQt version: {qtpy.PYQT_VERSION}")
    msgs.info(f"QtPy API_NAME: {qtpy.API_NAME}")

    defaultFont = app.font()
    msgs.info(f"Default font pixel size: {defaultFont.pixelSize()}")
    msgs.info(f"Default font point size: {defaultFont.pointSizeF()}")
    if defaultFont.pointSizeF() < 18.0:
        msgs.info(f"Setting font to 18.")
        defaultFont.setPointSize(18)
        app.setFont(defaultFont)


    main_window = MainWindow()
    main_window.setWindowTitle(main_window.tr("PypeIt Dashboard"))
    main_window.resize(1650,900)
    main_window.show()

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
