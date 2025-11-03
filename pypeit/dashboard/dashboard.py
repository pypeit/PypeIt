import sys
import subprocess
from pathlib import Path
from pypeit import msgs
import zmq

from qtpy.QtCore import QTimer, QSize, Qt, QMargins
from qtpy.QtGui import QIcon, QColor, QColorConstants, QPainter
from qtpy.QtWidgets import QApplication, QWidget, QHBoxLayout, QVBoxLayout, QPushButton, QGridLayout, QLabel, QProgressBar,QTabWidget, QListWidget, QAbstractItemView
import qtpy

# from pypeit.setup_gui.controller import start_gui
# from pypeit.scripts import setup

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
                           (QIcon.ThemeIcon.InputKeyboard, "Edit Setup", self.edit_setup_button),
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

        l = QLabel(text="keck_deimos_830g_m_8500.pypeit")
        l.setContentsMargins(value_cm)
        l.setStyleSheet(value_style_sheet)
        layout.addWidget(l,0,1,1,1)
        
        # --------------------- Calibration id group --------------
        calibration_id_label = QLabel(text="Calibration ID")
        layout.addWidget(calibration_id_label,1,0,1,1)

        l = QLabel(text="0")
        l.setContentsMargins(value_cm)
        l.setStyleSheet(value_style_sheet)
        layout.addWidget(l,1,1,1,1)#,alignment=Qt.AlignmentFlag.AlignLeft)
        
        # --------------------- Detector group ---------------------
        detector_label = QLabel(text="Detector")
        layout.addWidget(detector_label,2,0,1,1)

        l = QLabel(text="3")
        l.setContentsMargins(value_cm)
        l.setStyleSheet(value_style_sheet)
        layout.addWidget(l,2,1,1,1)
        
        # ------------------- Science file group ------------------
        science_file_label = QLabel(text="Science File")
        layout.addWidget(science_file_label,0,2,1,1)

        l = QLabel(text="DE.20100913.22358.fits.gz")
        l.setContentsMargins(value_cm)
        l.setStyleSheet(value_style_sheet)
        layout.addWidget(l,0,3,1,1)

        # ------------------- step group ----------------------
        step_label = QLabel(text="Step")
        layout.addWidget(step_label,1,2,1,1)

        l = QLabel(text="Calibrations")
        l.setContentsMargins(value_cm)
        l.setStyleSheet(value_style_sheet)
        layout.addWidget(l,1,3,1,1)

        # ------------------------ calibration step group ----------------
        calibration_step_label = QLabel(text="Calibration Step")
        layout.addWidget(calibration_step_label,2,2,1,1)

        l = QLabel(text="Tilts")
        l.setContentsMargins(value_cm)
        l.setStyleSheet(value_style_sheet)
        layout.addWidget(l,2,3,1,1)

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
        layout.addWidget(StatusWidget())
        tab_widget = QTabWidget()
        tab_widget.addTab(FilledBackgroundWidget(color=QColorConstants.Red),"QA")
        tab_widget.addTab(FileListWidget(),"Calibrations")
        tab_widget.addTab(FilledBackgroundWidget(color=QColorConstants.DarkGreen),"Science")

        layout.addWidget(tab_widget, 3)
        self.setLayout(layout)



# what does each button do. 
# the open setup button opens the setup at the current stage/ state
# edit setup allows you to enter a different setup file
# run all starts running the tasks and run does this step by step
# help does something

class MainWindow(QWidget):
    
    def __init__(self):
        super().__init__()

        layout = QHBoxLayout()
        setup_widget = ButtonWidget()
        dashboard_widget = DashboardWidget()
        layout.addWidget(setup_widget,alignment=Qt.AlignmentFlag.AlignTop)
        layout.addWidget(dashboard_widget,stretch=3)

        # -------- connections ---------
        setup_widget.open_setup_button.clicked.connect(self.start_controller)
        setup_widget.edit_setup_button.clicked.connect(self.edit_setup_file)

        self.setLayout(layout)
    def start_controller(self):
        command = ["pypeit_setup","--gui"]
        subprocess.Popen(command) # starting the controller runnner file
        # need to find the way to run this but it will be via a subprocess

        
    def edit_setup_file(self):
        file_path = b"setup_file=./sample_pypeit_file.pypeit"

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

    # subprocess.Popen([sys.executable, "-m", "controller_runner"]) # starting the controller runnner file
    # main_window.start_zmq()
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
