from functools import partial
import enum
from pathlib import Path

from qtpy.QtWidgets import QGroupBox, QHBoxLayout, QVBoxLayout, QComboBox, QToolButton, QFileDialog, QWidget, QGridLayout
from qtpy.QtWidgets import QMessageBox, QTabWidget, QTreeView, QTreeWidget, QTreeWidgetItem, QTableView, QPushButton, QProgressDialog, QDialog, QHeaderView
from qtpy.QtGui import QIcon
from qtpy.QtCore import Qt, Signal,QSettings, QStringListModel, QAbstractItemModel, QModelIndex

from pypeit.setup_gui.model import ModelState
from pypeit import msgs

class LocationPanel(QGroupBox):

    locationChanged = Signal(str)

    def __init__(self, parent, name, placeholder = None):
        super().__init__(title = name, parent=parent)
        self.parent = parent
        self.name = name
        self._settings = QSettings()
        self._settings.beginGroup(name.replace(" ", ""))
        history_list = self._settings.value("History")
        msgs.info(f"History: {repr(history_list)}, status: {repr(self._settings.status())}")
        self._history = QStringListModel()        
        if history_list is not None:
            if not isinstance(history_list, list):
                history_list = [history_list]
            self._history.setStringList(history_list)

        self._current_location = None
        layout = QHBoxLayout()
        layout.setSpacing(0)
        self.location = QComboBox(self)
        self.location.setEditable(True)
        self.location.setInsertPolicy(QComboBox.InsertAtTop)
        self.location.setModel(self._history)
        if placeholder is not None:
            self.location.lineEdit().setPlaceholderText(placeholder)

        self.location.setCurrentIndex(-1)
        self.location.activated.connect(self._new_location)
        layout.addWidget(self.location)
        self._browse_button = QToolButton()
        self._browse_button.setText("Browse...")
        self._browse_button.setIcon(QIcon.fromTheme("document-open"))
        self._browse_button.clicked.connect(self.browse)
        layout.addWidget(self._browse_button)
        self.setLayout(layout)

    @property
    def current_location(self):
        return self._current_location

    def browse(self):
        """Opens up a FileDialog requesting the user select a location.

        Returns: The new location, or None if Cancel was pressed.

        Signals: location.activated is emitted with the new location if 
        a new location is selected.
        """
        browse_dialog = QFileDialog(self.parent, caption = f"Select {self.name}")
        browse_dialog.setFileMode(QFileDialog.Directory)
        browse_dialog.setHistory(self._history.stringList())
        if browse_dialog.exec():
            selected_dir = browse_dialog.selectedFiles()[0]
            if selected_dir in self._history.stringList():
                self.location.setCurrentIndex(self._history.stringList().index(selected_dir))
            else:
                self.location.insertItem(0, selected_dir)
                self.location.setCurrentIndex(0)

            if selected_dir != self._current_location:
                self._current_location = selected_dir
                self.location.activated.emit(0)
            return self._current_location
        else:
            return None

    def _new_location(self, new_index):
        msgs.info(f"_new_location {new_index}")
        self.location.setCurrentIndex(new_index)        
        self._settings.setValue("History", self._history.stringList())
        self._current_location = self.location.currentText()
        #msgs.info(f"History: {repr(self._history.stringList())}")

    def setEnabled(self, value):
        self.location.setEnabled(value)
        self._browse_button.setEnabled(value)
        super().setEnabled(value)

    def updateLocation(self, new_value):
        msgs.info(f"location changed: {new_value}")
        if new_value is None and self._current_location is None or new_value == self._current_location:
            return
        self._current_location = new_value
        if new_value is None or new_value == "":
            self.location.setCurrentIndex(-1)
        elif new_value in self._history.stringList():
            self.location.setCurrentIndex(self._history.stringList().index(new_value))
        else:
            self.location.insertItem(0, new_value)
            self.location.setCurrentIndex(0)
            self._settings.setValue("History", self._history.stringList())
            #msgs.info(f"History: {repr(self._history.stringList())}")

class ObsLogTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)        
        group = QGroupBox("File Metadata")
        layout = QVBoxLayout()
        file_group_layout = QHBoxLayout()        
        self.obslog_table = PypeItSetupView.create_metadata_table(group)
        file_group_layout.addWidget(self.obslog_table)
        group.setLayout(file_group_layout)
        layout.addWidget(group)
        self.setLayout(layout)


    def setModel(self, model):
        self.obslog_table.setModel(model)
        self.obslog_table.horizontalHeader().setSortIndicator(model.sortColumn, model.sortOrder)
        self.obslog_table.setSortingEnabled(True)

class ConfigTab(QWidget):

    def __init__(self, model):
        super().__init__()

        layout = QVBoxLayout() 
        self._model = model
        file_group = QGroupBox("File Metadata")
        file_group_layout = QVBoxLayout()
        self.file_metadata_table = PypeItSetupView.create_metadata_table(file_group)
        self.file_metadata_table.setModel(model.metadata_model)
        file_group_layout.addWidget(self.file_metadata_table)
        file_group.setLayout(file_group_layout)        
        layout.addWidget(file_group)

        params_group = QGroupBox("PypeIt Parameters")
        params_group_layout = QHBoxLayout()
        self.params_tree = QTreeView(params_group)
        self.params_tree.setModel(model.params_model)
        self.params_tree.header().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.params_tree.expandToDepth(1)

        params_group_layout.addWidget(self.params_tree)
        params_group.setLayout(params_group_layout)        
        layout.addWidget(params_group)
        self.setLayout(layout)

    @property
    def name(self):
        return self._model.name

    @property
    def state(self):
        return self._model.state

class PypeItSetupView(QWidget):
    class Responses(enum.Enum):
        ACCEPT = enum.auto() # Also used for okkay, continue, etc
        SAVE   = enum.auto()
        CANCEL = enum.auto()


    def __init__(self):
        super().__init__(parent=None)

        self.layout = QGridLayout(self)

        self.layout.addWidget(self._create_options_panel(), 0, 0, 1, 1)
        self._setup_tab_panel = self._create_setup_tab_panel() 
        self.layout.addWidget(self._setup_tab_panel, 1, 0, 6, 1)
        self._setup_tab_panel.currentChanged.connect(self._currentTabChanged)
        self._create_button_box(7)

        # Todo get relative path to icon, also, is this really the icon we want to use?
        self.setWindowIcon(QIcon(str(Path(__file__).parent / "images/window_icon.png")))
        self.setWindowTitle("PypeIt Setup")

        self._config_tabs = []

    def display_error(self, message):
        QMessageBox.warning(self, "PypeIt Setup Error", message, QMessageBox.Ok)

    def prompt_for_save(self):
        response = QMessageBox.warning(self, "PypeIt Setup", "There are unsaved changes to PypeIt setup files.\nDo you want to save then?",
                                       QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel, QMessageBox.Save)
        if response == QMessageBox.Save:
            return PypeItSetupView.Responses.SAVE
        elif response == QMessageBox.Discard:
            return PypeItSetupView.Responses.ACCEPT
        else:
            return PypeItSetupView.Responses.CANCEL

    def start_operation_progress(self, op_name, max_progress_value, cancel_func):
        self.currentOpProgressDialog = QProgressDialog(op_name, "Cancel", 0, max_progress_value, parent=self)
        self.currentOpProgressDialog.setMinimumWidth(380)
        self.currentOpProgressDialog.setWindowTitle(op_name)
        self.currentOpProgressDialog.setMinimumDuration(1000)
        self.currentOpProgressDialog.setValue(0)            
        self.currentOpProgressDialog.canceled.connect(cancel_func)

    def set_operation_progress(self, increase, message=None):
        self.currentOpProgressDialog.setValue(self.currentOpProgressDialog.value() + increase)
        if message is not None:
            self.currentOpProgressDialog.setLabelText(message)

    def operation_complete(self):
        self.currentOpProgressDialog.done(QDialog.Accepted)
        self.currentOpProgressDialog = None

    def _create_options_panel(self):
        optionsPanel = QWidget(parent = self)
        layout = QHBoxLayout()

        spectrograph_box = QGroupBox(title="Spectrograph", parent=self)
        spectrograph_layout = QHBoxLayout()        

        self.spectrograph = QComboBox(spectrograph_box)
        self.spectrograph.setEditable(False)

        self.spectrograph.setPlaceholderText("select a spectrograph")
        spectrograph_layout.addWidget(self.spectrograph)
        spectrograph_box.setLayout(spectrograph_layout)
        layout.addWidget(spectrograph_box)

        self.raw_data = LocationPanel(self, "Raw Data Directory", placeholder="Choose raw data directory")
        layout.addWidget(self.raw_data)


        self.outdir = LocationPanel(self, "Output Directory", placeholder="Choose output directory")
        layout.addWidget(self.outdir)

        optionsPanel.setLayout(layout)
        return optionsPanel


    def _create_setup_tab_panel(self):
        tab_panel = QTabWidget(self)
        tab_panel.setTabPosition(QTabWidget.South)

        self.obslog = ObsLogTab(parent=tab_panel)
        tab_panel.addTab(self.obslog, "ObsLog")

        return tab_panel


    def setModel(self, model):
        self._model = model
        self.spectrograph.addItems(self._model.available_spectrographs)
        self.obslog.setModel(self._model.metadata_model)
        self._model.raw_data_dir_changed.connect(self.raw_data.updateLocation)
        self._model.configs_deleted.connect(self.deleteTabs)

    def createConfigTabs(self, config_models):
        self._setup_tab_panel.setUpdatesEnabled(False) # To prevent flickering when updating

        try:
            for config_model in config_models:
                # Keep tabs sorted
                tab_index = None
                for i in range(len(self._config_tabs)):
                    if config_model.name < self._config_tabs[i].name:
                        tab_index = i
                        break
                    elif config_model.name == self._config_tabs[i].name:
                        msgs.bug(f"Duplicate name in tab list: {config_model.name}")
                        return

                config_tab = ConfigTab(config_model)

                if config_model.state == ModelState.CHANGED:
                    tab_name = config_model.name + "*"
                else:
                    tab_name = config_model.name

                if tab_index == None:
                    tab_index = len(self._config_tabs)
                    self._config_tabs.append(config_tab)
                    self._setup_tab_panel.addTab(config_tab, tab_name)
                else:
                    self._config_tabs.insert(tab_index,config_tab)
                    # Insert at tab_index + 1 because of the ObsLog tab
                    self._setup_tab_panel.insertTab(tab_index + 1, config_tab, tab_name)

                config_model.stateChanged.connect(self._configStateChange)
        finally:
            self._setup_tab_panel.setUpdatesEnabled(True) # To allow redrawing after updating

    def get_current_tab_name(self):
        return self._setup_tab_panel.currentWidget().name

    def _currentTabChanged(self, index):
        if index == 0:
            # obslog
            self.saveTabButton.setEnabled(False)
        else:
            # Enable tab button if the model is in a changed state
            config_tab = self._config_tabs[index-1]
            self.saveTabButton.setEnabled(config_tab.state == ModelState.CHANGED)

    def _configStateChange(self, config_name):
        tab_index = [tab.name for tab in self._config_tabs].index(config_name)
        tab = self._config_tabs[tab_index]
        if tab.state == ModelState.CHANGED:
            self._setup_tab_panel.setTabText(tab_index+1, tab.name + "*")
        else:
            self._setup_tab_panel.setTabText(tab_index+1, tab.name)

    def deleteTabs(self, tab_list):
        msgs.info(f"View Deleting tabs {tab_list}")
        if len(tab_list) == 0:
            return

        # Go backwards through the config list so indices remain valid after removing
        for i in reversed(range(len(self._config_tabs))):
            if self._config_tabs[i].name in tab_list:
                # Remember to increment index by 1 because of the ObsLog tab
                self._setup_tab_panel.removeTab(i+1)
                del self._config_tabs[i]
    


    @classmethod
    def create_metadata_table(cls, parent):
        #default_columns =  ['filename', 'frametype', 'ra', 'dec', 'target', 'dispname', 'decker', 'binning', 'mjd', 'airmass', 'exptime']
        table = QTableView(parent)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        return table

    def _create_button_box(self, row):

        self.openButton = QPushButton("Open")
        self.clearButton = QPushButton("Clear")
        self.setupButton = QPushButton("Re-Run Setup")
        self.saveTabButton = QPushButton("Save Tab")
        self.saveAllButton = QPushButton("Save All")
        self.exitButton = QPushButton("Exit")

        buttonLayout = QHBoxLayout()
        self.layout.addLayout(buttonLayout, row, 0, 1, 1)
        buttonLayout.addWidget(self.openButton)
        buttonLayout.addWidget(self.clearButton)
        buttonLayout.addWidget(self.setupButton)
        buttonLayout.addWidget(self.saveTabButton)
        buttonLayout.addWidget(self.saveAllButton)
        buttonLayout.addStretch()
        buttonLayout.addWidget(self.exitButton)

    def promptForFile(self, caption, filter, history=[]):
        """Opens a dialog to prompt the user for an existing file.
        
        Args:
        caption (str): A caption for the dialog
        filter (str):  A QFileDialog filter for a file to open. For example: "Text Files (*.txt)"
        history (str): A list of file paths in the history of selected files. If this is not empty,
                       this list is given to the user as a list of previous files, and the last file
                       in the list is used as the starting directory for the dialog.

        Returns: The selected file, or None if Cancel was pressed.
        """
        file_dialog = QFileDialog(parent=self, caption = caption, filter=filter)        
        file_dialog.setFileMode(QFileDialog.ExistingFile)
        
        if len(history) > 0:
            file_dialog.setDirectory(history[-1])
            file_dialog.setHistory(history)

        if file_dialog.exec():
            selected_file = file_dialog.selectedFiles()[0]
            return selected_file
        else:
            return None

