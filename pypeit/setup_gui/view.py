from functools import partial
import enum
from pathlib import Path

from qtpy.QtWidgets import QGroupBox, QHBoxLayout, QVBoxLayout, QComboBox, QToolButton, QFileDialog, QWidget, QGridLayout
from qtpy.QtWidgets import QMessageBox, QTabWidget, QTreeView, QTreeWidget, QTreeWidgetItem, QTableView, QPushButton, QProgressDialog, QDialog, QHeaderView
from qtpy.QtGui import QIcon
from qtpy.QtCore import Qt, Signal,QSettings, QStringListModel, QAbstractItemModel, QModelIndex

from pypeit.setup_gui.model import ModelState
from pypeit import msgs


class PersistentStringListModel(QStringListModel):
    """
    Child class of QStringListModel that persists it's contents to QSettings.
    
    Args:
    settings_group (str):      The name of the settings group to save the string list under.
    settings_value_name (str): The name of the value to use when saving the string list.
    """
    def __init__(self, settings_group, settings_value_name):
        self._settings = QSettings()
        self._value_name = settings_value_name
        
        # Get persisted list (if any) and call superclass constructor with it
        self._settings.beginGroup(settings_group)
        list_value = self._settings.value(self._value_name)
        if list_value is None:
            list_value = []
        elif not isinstance(list_value, list):
            list_value = [list_value]

        super().__init__(list_value)

        # Attach to our own dataChanged signal to save when the list changes
        self.dataChanged.connect(self.save)

    def save(self, topLeft = None, bottomRight=None, roles=None):
        """
        saves the data persisted in this PersistentStringListModel. 

        Args:
        topLeft (QModelIndex):     The top left index of the data changhed.
                                   This is not and defaults to None.
        bottomRight (QModelIndex): The bottom right index of the data changed.
                                   This is not and defaults to None.
        roles (list):              List of roles changed. This is not and defaults to None.
        """
        # This receives the arguments from the dataChanged() signal 
        # but ignores them, persisting the entire list instead
        self._settings.setValue(self._value_name, self.stringList())

class LocationPanel(QGroupBox):
    """
    A custon widget that displays an editable combo box for file locations and a browse button
    to use a file dialog to enter the file location.  The history of past locations is kept in the 
    combo box list. 
    Args:
    parent(QWidget):      The parent object of the location panel.
    name (str):           The name to display for this location panel.
    browse_caption (str): The caption text to use when searching for locations, also used as place holder
                          text when no location has been set.

    Signals:
    location_changed(str): Signals when the location is changed.
    """

    location_changed = Signal(str)

    def __init__(self, parent=None, name=None, browse_caption = None):
        super().__init__(title = name, parent=parent)
        self.parent = parent
        self.name = name
        self.browse_caption = browse_caption
        # Setup the combo box
        layout = QHBoxLayout()
        layout.setSpacing(0)
        self._location = QComboBox(parent=self)
        self._location.setEditable(True)
        self._location.setInsertPolicy(QComboBox.InsertAtTop)

        # Setup history
        self._history = QStringListModel()
        self._location.setModel(self._history)
        if browse_caption is not None:
            self._location.lineEdit().setPlaceholderText(browse_caption)

        self._location.setCurrentIndex(-1)
        self._location.activated.connect(self._new_location)
        layout.addWidget(self._location)

        # Setup the browse button
        self._browse_button = QToolButton(parent=self)
        self._browse_button.setText(self.tr("Browse..."))
        self._browse_button.clicked.connect(self.browse)
        layout.addWidget(self._browse_button)
        self.setLayout(layout)

    @property
    def current_location(self):
        """The current location, or None if no location has been chosen."""
        return None if self._location.currentIndex()==-1 else self._location.currentText()

    def set_current_location(self, new_location):
        """Change the current location.
        
        Args:
        new_location (str): A new location. If this location is not already in the history, it will be added.
        """
        if new_location == self._location.currentText():
            # No change
            return

        elif new_location is None or new_location == "":
            # Set to place holder text
            self._location.setCurrentIndex(-1)
        else:
            # Set current text only sets the edit text, so we add it to the model if needed
            # and set the index
            idx = self._location.findText(new_location)
            if idx == -1:
                self._history.insertRows(0, 1, new_location)
                idx = 0
            self._location.setCurrentIndex(idx)

    def setModel(self, model):
        self._history = model
        self._location.setModel(model)
        self._location.setCurrentIndex(-1)

    def model(self):
        return self._history

    def browse(self, message=None):
        """Opens up a FileDialog requesting the user select a location.

        Returns: The new location, or None if Cancel was pressed.

        Signals: location.activated is emitted with the new location if 
        a new location is selected.
        """
        old_location = self.current_location
        browse_dialog = QFileDialog(self.parent, caption = self.name)
        browse_dialog.setFileMode(QFileDialog.Directory)
        browse_dialog.setHistory(self._history.stringList())
        if browse_dialog.exec():
            # the User selected a directory
            selected_dir = browse_dialog.selectedFiles()[0]

            if selected_dir in self._history.stringList():
                # If it's in the history, just select that entry (and don't insert a duplciate)
                self._location.setCurrentIndex(self._history.stringList().index(selected_dir))
            else:
                # It's not in the history, insert it at the top and select it
                self._location.insertItem(0, selected_dir)
                self._location.setCurrentIndex(0)

            # Use the new value as the current location, but don't
            # send a signal if it hasn't actually changed
            if selected_dir != old_location:
                self.location_changed.emit(selected_dir)
            return selected_dir
        else:
            return None

    def _new_location(self, new_index):
        """
        Signal handler for when the combo box selects a new location.

        Args:
        new_index: The index within the combo box that was selected.
        """
        # Forward the signal to clients with the actual string value of the new
        # location
        msgs.info(f"_new_location {new_index}")
        self.location_changed.emit(self.current_location)

    def setEnabled(self, value):
        """
        Set whether the widget (both combobox and browse button) is enabled.

        Args:
        value (bool): True to enable the widget is enabled, False to disable it
        """
        # This will also enable/disable the child combo box and button widgets
        super().setEnabled(value)

def _create_metadata_table(parent):
    """Method to create a table view for viewing file metadata.

    Args:
    parent (QtWidget): The parent widget of the table view
    """
    table = QTableView(parent)
    table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
    return table


class ObsLogTab(QWidget):
    """Widget for displaying the file metadata for all files in all configurations known to the GUI, aka the observation log.
    It is displayed as a tab within a QTabWidget within the PypeItSetupView object.

    Args:
    model (model.PypeItMetadataProxy): The model to use for the file metadta. If none, an
                                       empty tab is created.
    parent (QWidget):                  The parent widget of the tab.
    """

    def __init__(self, model=None, parent=None):
        super().__init__(parent)        
        group = QGroupBox(self.tr("File Metadata"))
        layout = QVBoxLayout()
        file_group_layout = QHBoxLayout()        
        self.obslog_table = _create_metadata_table(group)
        self.set_model(model)
        self.obslog_table.setSortingEnabled(True)

        file_group_layout.addWidget(self.obslog_table)
        group.setLayout(file_group_layout)
        layout.addWidget(group)
        self.setLayout(layout)

    def set_model(self, model):
        """
        Set the model represenging file metadata for the files in the obslog.
        
        Args:
        model (model.PypeItMetadataProxy): The model to use for file metadata.
        """
        if model is not None:
            self.obslog_table.setModel(model)
            self.obslog_table.horizontalHeader().setSortIndicator(model.sortColumn, model.sortOrder)


    @property
    def state(self):
        """
        Returns the state of this tab, always UNCHANGED  since we don't save the obs log
        """
        return ModelState.UNCHANGED

    @property
    def tab_name(self):
        """
        The name of this tab.

        Return (str): The name to use for this tab.
        """
        return "ObsLog"

class ConfigTab(QWidget):
    """Widget for displaying/editing the parameters and  file metadata for a PypeIt configuration. This is the information needed to 
    create a single .pypeit file. It is displayed as a tab within a QTabWidget within the PypeItSetupView object.

    Args: (model.SetupConfigurationModel): The model representing this configuration.
    """

    def __init__(self, model):
        super().__init__()

        layout = QVBoxLayout() 
        self._model = model

        # Create a group box and table view for the file metadata table
        file_group = QGroupBox(self.tr("File Metadata"))
        file_group_layout = QVBoxLayout()
        self.file_metadata_table = _create_metadata_table(file_group)
        self.file_metadata_table.setModel(model.metadata_model)
        file_group_layout.addWidget(self.file_metadata_table)
        file_group.setLayout(file_group_layout)        
        layout.addWidget(file_group)

        # Create a group box and a tree view for the pypeit parameters
        params_group = QGroupBox(self.tr("PypeIt Parameters"))
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
        """
        Return the name of this configuration.

        Return (str): The configuration name.
        """
        return self._model.name

    @property
    def state(self):
        """
        Return the state of this configuration.

        Return (model.ModelState): The state of this configuration's model. NEW, CHANGED, or UNCHANGED.
        """
        return self._model.state

    @property
    def tab_name(self):
        """
        The name of this tab. This may have a * in front of it if the tab has unsaved data.

        Return (str): The name to use for this tab.
        """
        if self.state == ModelState.CHANGED:
            return "*" + self.name
        else:
            return self.name

class PypeItSetupView(QTabWidget):
    """
    Tab widget which serves of the view of a PypeItSetup object.

    Args:
    model (model.PypeItSetupProxy): Proxy object adapting a PypeItSetup object
                                    to Qt.
    """
    def __init__(self, model):
        super().__init__()
        self._config_tabs = []
        self._obs_log_tab = ObsLogTab()
        self.addTab(self._obs_log_tab, self._obs_log_tab.tab_name)
        self._model = None
        self.setModel(model)
        self.setTabPosition(QTabWidget.South)

    def setModel(self, model):
        """
        Sets the model represenging a PypeItSetup object.

        Args:
        model (model.PypeItSetupProxy): The PypeItSetup model.
        """
        if self._model is not None:
            self._config_tabs = []
            self.clear()

        self._model = model
        self._obs_log_tab.set_model(model.metadata_model)
        self._model.configs_added.connect(self.create_config_tabs)
        self._model.configs_deleted.connect(self.delete_tabs)
        self.create_config_tabs(self._model.configs.values())

        
    def get_tab(self, index):
        """
        Return the tab at the given index.

        Args:
        index (int): Index of the tab to return.

        Return (ObsLogTab or ConfigTab): The tab at that index.
        """
        return self._obs_log_tab if index == 0 else self._config_tabs[index-1]

    @property
    def state(self):
        return ModelState.NEW if self._model is None else self._model.state

    def create_config_tabs(self, config_models):
        """
        Create new tabs for PypeIt configurations found in either a .pypeit file or
        a raw data directory.

        Args:
        config_models (list of model.SetupConfigurationModel): The new configuraitons.
        """
        msgs.info(f"create_config_tabs for {len(config_models)} configs")
        for config_model in config_models:
            try:
                msgs.info(f"Creating tab {config_model.name}")
                self.setUpdatesEnabled(False) # To prevent flickering when updating
                # Keep tabs sorted
                tab_index = None
                for i in range(len(self._config_tabs)):
                    if config_model.name < self._config_tabs[i].name:
                        tab_index = i
                        break
                    elif config_model.name == self._config_tabs[i].name:
                        msgs.bug(f"Duplicate name in tab list: {config_model.name}")
                        return

                if tab_index == None:
                    tab_index = len(self._config_tabs)

                config_tab = ConfigTab(config_model)
                config_model.stateChanged.connect(self.update_tab_status)

                self._config_tabs.insert(tab_index,config_tab)
                # Insert at tab_index + 1 because of the ObsLog tab
                self.insertTab(tab_index + 1, config_tab, config_tab.tab_name)

            finally:
                self.setUpdatesEnabled(True) # To allow redrawing after updating

    def delete_tabs(self, tab_list):
        """
        Delete any tabs no longer in the model.

        Args:
        tab_list (list of str): List of the configuration names removed.
        """
        msgs.info(f"View Deleting tabs {tab_list}")
        if len(tab_list) == 0:
            return

        # Go backwards through the config list so indices remain valid after removing
        for i in reversed(range(len(self._config_tabs))):
            if self._config_tabs[i].name in tab_list:
                # Remember to increment index by 1 because of the ObsLog tab
                self.removeTab(i+1)
                del self._config_tabs[i]
    
    def update_tab_status(self, name):
        """
        Update a tab's name to indicate if it has been changed but not saved.
        Args:
        tab_index (int): The index of the tab.
        """
        # Find the associated tab and set its name.
        for tab_index in range(len(self._config_tabs)):
            if self._config_tabs[tab_index].name == name:
                self.setTabText(tab_index+1, self._config_tabs[tab_index].tab_name)
                break

class DialogResponses(enum.Enum):
    """ Enum for the responses from a save dialog.
    Values:
    ACCEPT:  The dialog was accepted by the user. This could have been from an "Accept", "Ok", "Continue" or similar button.
    SAVE:    The user requested to save any unchanged data before continuing with an operation.
    CANCEL:  The user canceled the dialog."""
    ACCEPT = enum.auto()
    SAVE   = enum.auto()
    CANCEL = enum.auto()

class SetupGUIMainWindow(QWidget):
    """Main window widget for the PypeIt Setup GUI

    Args:
    available_spectrographs (list of str): The list of spectrographs supported by PypeIt.
    controller (SetupGUIController): The controller for the PypeITSetupGUI.
    """


    def __init__(self, pypeit_setup, controller):
        super().__init__(parent=None)

        self.layout = QGridLayout(self)
        self.pypeit_setup = pypeit_setup    
        self.controller = controller

        # Layout the main window
        self.layout.addWidget(self._create_options_panel(), 0, 0, 1, 1)
        self.setup_view = PypeItSetupView(pypeit_setup) 
        self.setup_view.currentChanged.connect(self.current_tab_changed)
        self.layout.addWidget(self.setup_view, 1, 0, 6, 1)
        self._create_button_box(7)

        # Setup application/window icon TODO this doesn't work in windows. Mac???
        self.setWindowIcon(QIcon(str(Path(__file__).parent / "images/window_icon.png")))
        self.setWindowTitle(self.tr("PypeIt Setup"))

        self.resize(1650,900)

    def display_error(self, message):
        """
        Display an error message pop up dialog to the user.

        Args:
        message (str): The message to display.
        """
        QMessageBox.warning(self, "PypeIt Setup Error", message, QMessageBox.Ok)

    def prompt_for_save(self):
        """
        Prompt the user if they want to save any unsaved data.

        Return (view.DialogResponses):  
            SAVE if the user wants to save, ACCEPT if they want to continue without 
            saving (i.e. discard the unsaved data), CANCEL if they want to cancel the current
            operation and to not lose any data.
        """
        response = QMessageBox.warning(self, self.tr("PypeIt Setup"), self.tr("There are unsaved changes to PypeIt setup files.\nDo you want to save then?"),
                                       QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel, QMessageBox.Save)
        if response == QMessageBox.Save:
            return DialogResponses.SAVE
        elif response == QMessageBox.Discard:
            return DialogResponses.ACCEPT
        else:
            return DialogResponses.CANCEL

    def start_operation_progress(self, op_caption, max_progress_value, cancel_func):
        """Start displaying progress information for an operation. This uses the QProgressDialog, which will not
        display itself until a minimum amount of time has passed (currently 1s)
        
        Args:
        op_caption (str):         The name of the operation.
        max_progress_value (int): The maximum progress value (i.e. the value when done).
        cancel_func (callable):   A callable to deal with cancel being pressed in the 
                                  progress dialog.
        """
        self.current_op_progress_dialog = QProgressDialog(self.tr(op_caption), self.tr("Cancel"), 0, max_progress_value, parent=self)
        self.current_op_progress_dialog.setMinimumWidth(380)
        self.current_op_progress_dialog.setWindowTitle(op_caption)
        self.current_op_progress_dialog.setMinimumDuration(1000)
        self.current_op_progress_dialog.setValue(0)            
        self.current_op_progress_dialog.canceled.connect(cancel_func)

    def increase_operation_progress(self, increase, message=None):
        """
        Increase the amount of progress being displayed in a progress dialog)
        Args:
        incrase (int): How much to increase the current progress by.
        message (str): A message indicating what step has been performed. (Optional)
        """
        self.current_op_progress_dialog.setValue(self.current_op_progress_dialog.value() + increase)
        if message is not None:
            self.current_op_progress_dialog.setLabelText(message)

    def operation_complete(self):
        """
        Stop displaying progress for an operation because it has completed..
        """
        self.current_op_progress_dialog.done(QDialog.Accepted)
        self.current_op_progress_dialog = None
        self.save_all_button.setEnabled(self.setup_view.state==ModelState.CHANGED)
        self.clearButton.setEnabled(self.setup_view.state!=ModelState.NEW)

    def _create_options_panel(self):
        """
        Create the panel with the Raw Data Directory, spectrograph, and output file
        directory options.

        Returns (QtWidget): The new panel.
        """
        optionsPanel = QWidget(parent = self)
        layout = QHBoxLayout()

        spectrograph_box = QGroupBox(title=self.tr("Spectrograph"), parent=self)
        spectrograph_layout = QHBoxLayout()        

        self.spectrograph = QComboBox(spectrograph_box)
        self.spectrograph.setEditable(False)
        self.spectrograph.addItems(self.pypeit_setup.available_spectrographs)
        self.spectrograph.setPlaceholderText("select a spectrograph")
        self.spectrograph.setCurrentIndex(-1)
        self.spectrograph.currentIndexChanged.connect(self.current_spec_changed)
        self.pypeit_setup.spectrograph_changed.connect(self.set_spectrograph)

        spectrograph_layout.addWidget(self.spectrograph)
        spectrograph_box.setLayout(spectrograph_layout)
        layout.addWidget(spectrograph_box)


        # Create Raw Data Directory location panel with persistent history
        self.raw_data = LocationPanel(self, self.tr("Raw Data Directory"), browse_caption=self.tr("Choose raw data directory"))
        self.raw_data.setModel(PersistentStringListModel("RawDataDirectory", "History"))
        self.raw_data.setEnabled(False)
        self.raw_data.location_changed.connect(self.run_setup)
        self.pypeit_setup.raw_data_dir_changed.connect(self.set_raw_data_dir)
        layout.addWidget(self.raw_data)

        # Create Output Directory location panel with persistent history
        self.outdir = LocationPanel(self, self.tr("Output Directory"), browse_caption=self.tr("Choose output directory"))
        self.outdir.setModel(PersistentStringListModel("OutputDirectory", "History"))
        layout.addWidget(self.outdir)

        optionsPanel.setLayout(layout)
        return optionsPanel            
        
    def current_tab_changed(self, index):
        """
        Handles enabling/disabling the save tab button when the current tab in the
        setup view tab panel changes.

        Args:
        index (int): The
        """
        # Enable tab button if the model is in a changed state
        config_tab = self.setup_view.get_tab(index)

        # Enable tab button if the model is in a changed state
        self.saveTabButton.setEnabled(config_tab.state == ModelState.CHANGED)

    def set_spectrograph(self, new_spectrograph):
        """
        Updates the spectrograph combobox to the new value. 

        Args:
        new_spectrograph (str): The new spectrograph value. It's the caller's responsibility to validate that this is
                                a valid spectrograph. If it is not, the combo box will reset to empty.
        """
        self.spectrograph.setCurrentText(new_spectrograph)
        if self.spectrograph.currentIndex() == -1:
            # No spectrograph set
            self.raw_data.setEnabled(False)
            self.setupButton.setEnabled(False)
        else:
            self.raw_data.setEnabled(True)

    def set_raw_data_dir(self, new_raw_data_dir):
        """
        Updates the raw data directory LocationPanel to the new value. 

        Args:
        new_spectrograph (str): The new location.
        """
        self.raw_data.set_current_location(new_raw_data_dir)
        if self.raw_data.current_location is None:
            self.setupButton.setEnabled(False)
        else:
            self.setupButton.setEnabled(True)

    def current_spec_changed(self):
        """
        Handles enabling/disabling the Raw Data Directory combo box when the spectrograph is set, 
        notifies controller of new spectrograph.
        """
        self.controller.spectrograph_changed(self.spectrograph.currentText())

    def run_setup(self):
        """
        Starts running setup on the currently selected spectrograph and raw data directory.
        """
        self.controller.run_setup(self.raw_data.current_location)

    def save_tab(self):
        """
        Saves the currently selected tab.
        """
        current_tab_idx = self.setup_view.currentIndex()
        # Don't save if there are no tabs (empty PypeItSetup or the obs log is the current tab)
        if current_tab_idx > 0:
            # Make sure there's an output location
            location = self.outdir.current_location
            if location is None or location == "":
                location=self.outdir.browse()
                if location is None:
                    # Prompt was canceled
                    return
            tab = self.setup_view.get_tab(current_tab_idx)
            self.controller.save_tab(tab.name, self.outdir.current_location)
            self.saveTabButton.setEnabled(tab.state == ModelState.CHANGED)
            self.save_all_button.setEnabled(self.setup_view.state==ModelState.CHANGED)

    def save_all(self):
        """
        Saves all unsaved tabs.
        """
        # Make sure there's an output location
        location = self.outdir.current_location
        if location is None or location == "":
            location=self.outdir.browse()
            if location is None:
                # Prompt was canceled
                return

        self.controller.save_all(location)
        self.save_all_button.setEnabled(self.setup_view.state==ModelState.CHANGED)

    def clear(self):
        self.controller.clear()
        self.save_all_button.setEnabled(self.setup_view.state==ModelState.CHANGED)
        self.clearButton.setEnabled(self.setup_view.state!=ModelState.NEW)

    def open_pypeit_file(self):
        self.controller.open_pypeit_file()
        self.save_all_button.setEnabled(self.setup_view.state==ModelState.CHANGED)
        self.clearButton.setEnabled(self.setup_view.state!=ModelState.NEW)


    def _create_button_box(self, row):
        """Create the box with action buttons.

        Args:
        row (int): The row to add the button box to within the
                   main windows layout.
        callback (callable): A callback to receive notification of when a button
                             starts an operation. This callback should take one keyword
                             value named "op" of type controller.Operation.

        """
            
        button_layout = QHBoxLayout()
        self.layout.addLayout(button_layout, row, 0, 1, 1)

        button = QPushButton(text = 'Open')
        button.setToolTip("Open a .pypeit file.")
        button.clicked.connect(self.open_pypeit_file)
        button_layout.addWidget(button)

        button = QPushButton(text = 'Clear')
        button.setToolTip("Clear everything and start with a blank slate.")
        button.setEnabled(False)
        button.clicked.connect(self.clear)
        button_layout.addWidget(button)
        self.clearButton = button


        button = QPushButton(text = 'Re-Run Setup')
        button.setToolTip("Scan the raw data and setup a .pypeit file for unique configuration.")
        button.setEnabled(False)
        button.clicked.connect(self.run_setup)
        button_layout.addWidget(button)
        self.setupButton = button

        button = QPushButton(text = 'Save Tab')
        button.setToolTip("Save the curretly active tab.")
        button.setEnabled(False)
        button.clicked.connect(self.save_tab)
        button_layout.addWidget(button)
        self.saveTabButton = button

        button = QPushButton(text = 'Save All')
        button.setToolTip("Save all tabs with unsaved changes.")
        button.setEnabled(False)
        button.clicked.connect(self.save_all)
        button_layout.addWidget(button)
        self.save_all_button = button

        button_layout.addStretch()

        button = QPushButton(text = 'Exit')
        button.setToolTip("Quits this application.")
        button.clicked.connect(self.controller.exit)
        button_layout.addWidget(button)

    def prompt_for_file(self, caption, filter, history=[]):
        """Opens a dialog to prompt the user for an existing file.
        
        Args:
        caption (str): A caption for the dialog
        filter (str):  A QFileDialog filter for a file to open. For example: "Text Files (*.txt)"
        history (str): A list of file paths in the history of selected files. If this is not empty,
                       this list is given to the user as a list of previous files, and the last file
                       in the list is used as the starting directory for the dialog.

        Returns: The selected file, or None if Cancel was pressed.
        """
        file_dialog = QFileDialog(parent=self, caption = self.tr(caption), filter=filter)        
        file_dialog.setFileMode(QFileDialog.ExistingFile)
        
        if len(history) > 0:
            file_dialog.setDirectory(history[-1])
            file_dialog.setHistory(history)

        if file_dialog.exec():
            selected_file = file_dialog.selectedFiles()[0]
            return selected_file
        else:
            return None

