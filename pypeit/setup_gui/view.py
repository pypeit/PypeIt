from functools import partial
import enum
from pathlib import Path

from qtpy.QtWidgets import QGroupBox, QHBoxLayout, QVBoxLayout, QComboBox, QToolButton, QFileDialog, QWidget, QGridLayout, QFormLayout
from qtpy.QtWidgets import QMessageBox, QTabWidget, QTreeView, QLayout, QLabel, QScrollArea, QListWidget, QTableView, QPushButton, QProgressDialog, QDialog, QHeaderView, QSizePolicy, QCheckBox, QDialog
from qtpy.QtGui import QIcon, QPalette, QColor, QValidator
from qtpy.QtCore import Qt, QSize, Signal,QSettings, QStringListModel, QAbstractItemModel, QModelIndex, QMargins

from pypeit.setup_gui.model import ModelState, ObservingConfigModel, available_spectrographs
from pypeit import msgs


class DialogResponses(enum.Enum):
    """ Enum for the responses from a dialog."""

    ACCEPT         = enum.auto()
    """The user accepted the dialog. This could have been from an "Accept", "Ok", "Continue" or similar button."""

    ACCEPT_FOR_ALL = enum.auto()
    """The user accepted the dialog for this and any subsequent requests. For example, a Save All that should save everthing to the same path."""

    SAVE           = enum.auto()
    """The user requested to save any unchanged data before continuing with an operation."""

    CANCEL         = enum.auto()
    """The user canceled the dialog."""


class FileDialog():
    """Opens a file dialog for either opening or saving files.  This encapsulates
    the logic to use QFileDialog in a way that can be used by the view or controller and
    can be easilly patched by unit tests.
    
    Args:
        parent (QWidget):                  Parent of the dialog widget.
        caption (str):                     Caption for the dialog.
        file_mode (QFileDialog.FileMode):  The file mode for the dialog.
        filter (str):                      The filter to use for files in the file dialog. For example: "PypeIt input files (\*.pypeit)"
        history (list of str):             The list of paths to populate the history of the dialog.
        save_all (bool):                   Whether the dialog should present a "Use this location for everything" option
                                           when saving.
     """
    def __init__(self, parent, caption, file_mode, filter=None, history=[], save_all=False):
        self._dialog = QFileDialog(parent, caption=parent.tr(caption), filter=filter)
        self._save_all = save_all

        self._dialog.setOption(QFileDialog.Option.DontUseNativeDialog)
        self._dialog.setFileMode(file_mode)
        if len(history) > 0:
            self._dialog.setDirectory(history[-1])
            self._dialog.setHistory(history)

        # Adding a checkbox like this only works for non-native dialogs. And if the 
        # default layout of the QFileDialog ever changes this could stop working.
        if self._save_all:
            sd_layout = self._dialog.layout()
            row_count = sd_layout.rowCount()
            self.use_for_all_checkbox = QCheckBox(self._dialog.tr("Use this location for everthing."), self._dialog)
            sd_layout.addWidget(self.use_for_all_checkbox, row_count, 0)

        self.selected_path = ""

    def show(self):
        """Show a modal file dialog. If the dialog is accepted, the result will be saved
        in the selected_path attribute.
        
        Returns: 
            DialogResponses: CANCEL, ACCEPT_FOR_ALL, or ACCEPT.
        """
        result = self._dialog.exec()
        if result ==  QDialog.Rejected:
            return DialogResponses.CANCEL
        else:
            self.selected_path = self._dialog.selectedFiles()[0]
            if self._save_all and self.use_for_all_checkbox.isChecked():
                return DialogResponses.ACCEPT_FOR_ALL
            else:
                return DialogResponses.ACCEPT


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
        Saves the data persisted in this PersistentStringListModel. 

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
    A custon widget that displays an editable combo box for entering file locations, a browse button
    to use a file dialog to enter the file location, and a list of all the file locations entered.
    The history of past locations is kept in the combo box list and file dialog history. 

    Args:
        parent(QWidget):                 The parent object of the location panel.
        name (str):                      The name to display for this location panel.
        browse_caption (str):            The caption text to use when searching for locations, also used as place holder
                                         text when no location has been set.
        lines_to_display (int,Optional): How many lines to display in the list of file locations. Defaults to 5.
    """

    location_added = Signal(str)
    """Signal(str): Signals when a location is added."""

    def __init__(self, parent=None, name=None, browse_caption = None, lines_to_display=5):
        super().__init__(title = name, parent=parent)
        self.parent = parent
        self.name = name
        self.browse_caption = browse_caption
        # Setup the combo box
        layout = QGridLayout()
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
        layout.addWidget(self._location, 0, 0, 1, 1)

        # Setup the browse button
        self._browse_button = QToolButton(parent=self)
        self._browse_button.setText(self.tr("Browse..."))
        self._browse_button.clicked.connect(self.browse)
        layout.addWidget(self._browse_button, 0, 1, 1, 1)

        # The list of current locations
        self._location_list = QListWidget(self)        
        ll_margins = self._location_list.contentsMargins()
        ll_fm = self._location_list.fontMetrics()
        self._location_list.setFixedHeight(ll_margins.top() + ll_margins.bottom() + ll_fm.height()*lines_to_display)

        layout.addWidget(self._location_list, 1, 0, 1, 2)

        self.setLayout(layout)

    @property
    def locations(self):
        """The list of current locations."""
        return [self._location_list.item(i).text() for i in range(self._location_list.count())]

    def set_locations(self, new_locations):
        """Change the current list of locations.
        
        Args:
            new_locations (list of str): The list of new locations
        """
        self._location_list.clear()
        self._location_list.addItems(new_locations)

    def add_location(self, new_location):
        """Add a new location to the list.

        Args:
            new_location (str): The new location to add.
        """

        # Add to history if needed
        if new_location not in self._history.stringList():
            # It's not in the history, insert it into the 
            # location list (which uses the history as its model)
            self._location.insertItem(0, new_location)

        # Don't add duplicate items
        items = self._location_list.findItems(new_location, Qt.MatchExactly)
        if items is None or len(items) == 0:
            self._location_list.addItem(new_location)
            self.location_added.emit(new_location)

    def setHistory(self, history):
        self._history = history
        self._location.setModel(history)
        self._location.setCurrentIndex(-1)

    def history(self):
        return self._history

    def browse(self):
        """Opens up a :class:`FileDialog` requesting the user select a location.

        Returns: 
            str: The new location, or None if Cancel was pressed.
        """
        browse_dialog = FileDialog(self.parent, 
                                   caption = self.name, 
                                   file_mode=QFileDialog.Directory,
                                   history=self._history.stringList())

        if browse_dialog.show() == DialogResponses.ACCEPT:
            # the User selected a directory
            # Add it to our location list
            self.add_location(browse_dialog.selected_path)

            return browse_dialog.selected_path
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
        if new_index != -1: 
            new_location = self._location.currentText()
            self.add_location(new_location)
            self._location.setCurrentIndex(-1)

    def setEnabled(self, value):
        """
        Set whether the widget (both combobox and browse button) is enabled.

        Args:
            value (bool): True to enable the widget is enabled, False to disable it
        """
        # This will also enable/disable the child combo box and button widgets
        super().setEnabled(value)

class PypeItMetadataView(QTableView):
    """QTableView to display file metadata.

    Args:
        parent (QtWidget):                                            The parent widget of the table view
        model  (:class:`pypeit.setup_gui.model.PypeItMetadataProxy`): The model for the table. Optional, defaults to None.
    """
    def __init__(self, parent, model=None):
        super().__init__(parent=parent)
        self.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.setModel(model)
    
    def setModel(self, model):
        """Set the PypeItMetadataProxy model to use for the table.
        
        Args:
            model (:class:`pypeit.setup_gui.model.PypeItMetadataProxy`):  The model for the table, None for no model.
        """
        super().setModel(model)
        if model is None:
            self.setSortingEnabled(False)
        else:
            self.setSortingEnabled(True)
            self.horizontalHeader().setSortIndicator(model.sortColumn, model.sortOrder)


class ConfigValuesPanel(QGroupBox):

    """Scrollable panel to display configuration values for one or more frames.
    
    Args:
        spectrograph (str):         Name of spectrograph for the configuration.
        config (list of tuple):     The name/value pairs for the configuration keys defined by the spectrograph.
        lines_to_display (int):     How many lines to display before scrolling.
        parent (QWidget, Optional): The parent widget, defaults to None.
    """
    def __init__(self, spectrograph, config, lines_to_display, parent=None):

        super().__init__(parent=parent)        
        self.setTitle(self.tr("Setup"))

        # Put everything in a scroll area
        layout = QHBoxLayout(self)
        scroll_area = QScrollArea(parent=self)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        # A widget using a form layout to display the spectrograph + configuration values
        form_widget = QWidget()
        form_widget_layout = QFormLayout(form_widget)

        # Keep this from expanding too large.
        form_widget_layout.setSizeConstraint(QLayout.SizeConstraint.SetMinimumSize)
        
        # First line is the spectrograph.
        form_widget_layout.addRow(self.tr("Spectrograph"), QLabel(spectrograph))

        # Add additional rows for configuration keys
        for key, value in config:
            form_widget_layout.addRow(key, QLabel(str(value)))

        # Don't add extra margins in the FormLayout
        form_widget_layout.setContentsMargins(0, 0, 0, 0)
        scroll_area.setWidget(form_widget)

        # Figure out the correct height for this panel, so that only the spectrograph and self.number_of_lines
        # config keys are visible

        # Find the minimum height of the form widget needed to hold the number of lines to display
        fm = self.fontMetrics()
        min_fw_height = form_widget_layout.verticalSpacing()*(lines_to_display-1) + fm.height()*lines_to_display

        # The height of this panel is that height plus the margins + the group box title
        scroll_area_margins = scroll_area.contentsMargins()
        group_box_margins = self.contentsMargins()
        form_widget_margins = form_widget.contentsMargins()
        self.setFixedHeight(min_fw_height + 
                            fm.height()   +  # Group Box Title
                            group_box_margins.top()   + group_box_margins.bottom() +
                            scroll_area_margins.top() + scroll_area_margins.bottom() +
                            form_widget_margins.top() + form_widget_margins.bottom())

        # Set horizontal policy to use our minimum size, and vertical to use the fixed size we just set
        self.setSizePolicy(QSizePolicy(horizontal=QSizePolicy.Minimum, vertical = QSizePolicy.Fixed, type=QSizePolicy.DefaultType))
        
        layout.addWidget(scroll_area)

class PypeItFileTab(QWidget):
    """Widget for displaying the information needed for one pypeit file. This includes
    the spectrograph, the configuration keys and values, the files that share that 
    configuration, and the PypeIt parameters.

    Args:
        (pypeit.setup_gui.model.PypeItFileModel): The model representing all the information needed for a .pypeit file.
    """

    def __init__(self, model):
        super().__init__()

        layout = QVBoxLayout(self) 
        self._model = model

        # Add the file name as the first row and second rows
        filename_label = QLabel(self.tr("Filename:"))
        filename_label.setAlignment(Qt.AlignLeft)
        layout.addWidget(filename_label)

        self.filename_value = QLabel(model.filename)
        self.filename_value.setAlignment(Qt.AlignLeft)
        layout.addWidget(self.filename_value)

        # Add the spectrograph configuration keys to the third row, first column
        third_row_layout = QHBoxLayout()
        layout.addLayout(third_row_layout)

        # Build a list of the configuration keys, to avoid displaying them in the
        # arbitrary order chosen by the dict
        config = [(key, model.config.config_dict[key]) for key in model.config.get_config_keys()]

        # Add the ConfigValuesPanel, displaying the spectrograph + config keys.
        config_panel = ConfigValuesPanel(model.config.spectrograph, config, 5, parent=self)
        third_row_layout.addWidget(config_panel)

        # Add the Raw Data directory panel to the third row, second column
        self.raw_data_paths = LocationPanel(self, self.tr("Raw Data Directories"), browse_caption=self.tr("Choose raw data directory"), lines_to_display=5)
        self.raw_data_paths.setHistory(PersistentStringListModel("RawDataDirectory", "History"))
        self.raw_data_paths.set_locations(model.metadata_model.get_metadata_paths())

        third_row_layout.addWidget(self.raw_data_paths)

        # Make the metadata wider than the config panel
        third_row_layout.setStretch(1, 2)

        # Create a group box and table view for the file metadata table
        file_group = QGroupBox(self.tr("File Metadata"))
        file_group_layout = QVBoxLayout()
        self.file_metadata_table = PypeItMetadataView(self, model.metadata_model)
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
        
        # Stretch the metadata and params rows more than the filename and config_key rows
        layout.setStretch(2,4)
        layout.setStretch(3,10)
        layout.setStretch(4,10)

        self._model.state_changed.connect(self.update_filename)

    """
    Signal handler that updates the filename when the underlying model is saved.
    """
    def update_filename(self):
        self.filename_value.setText(self._model.filename)

    @property
    def name(self):
        """
        str: The configuration name.
        """
        return self._model.config.name

    @property
    def state(self):
        """
        :class:`pypeit.setup_gui.model.ModelState`): The state of this configuration's model. NEW, CHANGED, or UNCHANGED.
        """
        return self._model.state

    @property
    def tab_name(self):
        """
        str: The name to display for this tab. The name will  have a "*" in it if it has been modified.
        """
        if self.state == ModelState.CHANGED:
            return "*" + self.name
        else:
            return self.name

class ObsLogTab(QWidget):
    """Widget for displaying the observation log for raw data files for the same spectrograph but 
    potentially different observing configurations.

    Args:
        model (:class:`pypeit.setup_gui.model.PypeItSetupModel`): The model to use for the file metadata.
        parent (QWidget): The parent widget of the tab.
    """

    def __init__(self, model, parent=None):
        super().__init__(parent)        

        layout = QVBoxLayout(self)
        # Place the spectrograph group box and combo box in the first row, first column
        top_row_layout = QHBoxLayout()
        layout.addLayout(top_row_layout)
        spectrograph_box = QGroupBox(title=self.tr("Spectrograph"), parent=self)
        spectrograph_layout = QHBoxLayout()        

        self.spectrograph = QComboBox(spectrograph_box)
        self.spectrograph.addItems(available_spectrographs())
        self.spectrograph.setCurrentIndex(-1)
        self.spectrograph.setEditable(True)
        self.spectrograph.lineEdit().setPlaceholderText(self.tr("Select a spectrograph"))
        self.spectrograph.setInsertPolicy(QComboBox.NoInsert)
        self.spectrograph.setValidator(SpectrographValidator())
        spectrograph_layout.addWidget(self.spectrograph)
        spectrograph_layout.setAlignment(self.spectrograph, Qt.AlignTop)
        spectrograph_box.setLayout(spectrograph_layout)
        top_row_layout.addWidget(spectrograph_box)

        # Add the Raw Data directory panel to the first row, second column
        self.raw_data_paths = LocationPanel(self, self.tr("Raw Data Directories"), browse_caption=self.tr("Choose raw data directory"), lines_to_display=5)
        self.raw_data_paths.setHistory(PersistentStringListModel("RawDataDirectory", "History"))
        self.raw_data_paths.setEnabled(False)

        top_row_layout.addWidget(self.raw_data_paths)
        # Make the metadata wider than the spectrograph
        top_row_layout.setStretch(1, 2)


        # Add the File Metadata box in the second row
        file_group = QGroupBox(self.tr("File Metadata"))
        file_group_layout = QHBoxLayout()        
        self.obslog_table = PypeItMetadataView(file_group, model.metadata_model)
        file_group_layout.addWidget(self.obslog_table)
        file_group.setLayout(file_group_layout)
        layout.addWidget(file_group)
        # Make File Metadata taller than the spectrograph/raw data paths row
        layout.setStretch(1,4)
        self.setModel(model)

        # Update model with new spectrograph/data paths
        self.spectrograph.textActivated.connect(self._model.set_spectrograph)
        self.spectrograph.textActivated.connect(self.update_raw_data_paths_state)
        self.raw_data_paths.location_added.connect(self._model.add_raw_data_directory)

    @property
    def state(self):
        """pypeit.setup_gui.model.ModelState: The state of this tab, either "NEW" if there's no metadata or UNCHANGED  
        if there is.  It can never be "CHANGED" because we don't save the obs log.
        """
        return ModelState.UNCHANGED

    @property
    def tab_name(self):
        """str: The name of this tab, always "ObsLog". Included for consistency with :class:`PypeItFileTab`.
        """
        return "ObsLog"

    @property
    def name(self):
        """str: The name of this unique configuration.  Always "ObsLog" for the ObsLog tab.
        Included for consistency with :class:`PypeItFileTab`."""
        return "ObsLog"

    def setModel(self,model):
        """Set a new model for file metadata.

        Args:
            model (:class:`pypeit.setup_gui.model.PypeItSetupModel`): The new metadata model
        """
        self._model=model
        self.raw_data_paths.set_locations(model.metadata_model.get_metadata_paths())
        self.obslog_table.setModel(model.metadata_model)
        
        # Update based on changes to the model
        # We don't listen specifically to spectrograph/raw data dir changes because
        # those are sent when our widgets are updated, and we don't want an infinite loop
        # of signals
        model.state_changed.connect(self.update_from_model)

    def update_raw_data_paths_state(self):
        """Enable/Disable the raw data paths location panel based on the model state. """
        if self._model.state == ModelState.NEW:
            if self.spectrograph.currentIndex() == -1:
                self.raw_data_paths.setEnabled(False)
            else:
                self.raw_data_paths.setEnabled(True)
        else:
            self.raw_data_paths.setEnabled(True)

    def update_from_model(self):
        """
        Updates the spectrograph and raw data location widgets based on model updates.
        """
        # Update the current spectrograph
        if self._model.spectrograph is None:
            self.spectrograph.setCurrentIndex(-1)
        else:
            self.spectrograph.setCurrentText(self._model.spectrograph)

        # Disable changing the spectrograph if the model isn't in the new state
        if self._model.state == ModelState.NEW:
            self.spectrograph.setEnabled(True)
        else:
            self.spectrograph.setEnabled(False)

        # Get the latest raw data dirs
        self.raw_data_paths.set_locations(self._model.metadata_model.get_metadata_paths())
        self.update_raw_data_paths_state()


class SpectrographValidator(QValidator):
    """Validator to check whether a spectrograph name is valid, or potentially valid.
    This is used by the spectrograph combo box to allow tab completion without
    allowing invalid names to be entered."""
    def __init__(self):
        self._supported_spectrographs = available_spectrographs()
        super().__init__()

    def validate(self, str_input, int_input):
        """Validates whether a string is a valid spectrograph name, a
        prefix of a valid spectrograph name, or invalid. The comparison is
        case insensitive because the combo box QCompleter object will convert it to
        the appropriate case. Overriden method from QValidator.
        
        Args:
            str_input: A string input to validate.
            int_input: An integer input to validate (not used).

        Returns:
            QValidator.State: Acceptable, Intermediate, or Invalid based on the input.
        """
        if str_input is not None:
            if str_input.lower() in self._supported_spectrographs:
                return QValidator.Acceptable
            else:
                for spectrograph in self._supported_spectrographs:
                    if spectrograph.startswith(str_input.lower()):
                        return QValidator.Intermediate
        return QValidator.Invalid

class PypeItSetupView(QTabWidget):
    """
    Widget which serves of the view of a :class:`PypeItSetup` object.
    It allows the user to specify file paths and a spectrograph,
    and using tabs for each unique configuration that can be
    written to a PypeIt file.

    Args:
        model (:class:`pypeit.setup_gui.model.PypeItSetupModel`): Model object for a PypeItSetup.
    """

    def __init__(self, model):
        super().__init__()
        # Setup our tabs
        self._file_tabs = []
        self._obs_log_tab = ObsLogTab(model)

        # Add the tab widget for the obs log and pypeit files in the second row
        self.addTab(self._obs_log_tab, self._obs_log_tab.tab_name)
        self.setTabPosition(QTabWidget.South)


        # Set our model
        self._model = model
        self._obs_log_tab.setModel(self._model)

        self.create_file_tabs(self._model.pypeit_files.values())

        # Get notifications about new/removed unique configurations
        self._model.configs_added.connect(self.create_file_tabs)
        self._model.configs_deleted.connect(self.delete_tabs)

    @property
    def state(self):
        """:class:`pypeit.setup_gui.model.ModelState`: The state of the underlying model."""
        return ModelState.NEW if self._model is None else self._model.state

    def create_file_tabs(self, pypeit_file_models):
        """
        Create new tabs for new unique configurations found either by openeing a pypeit file or
        by reading raw data directories.

        Args:
            pypeit_file_models (list of :class:`pypeit.setup_gui.model.PypeItFileModel`): Models for the tabs to add.
        """
        msgs.info(f"create_file_tabs for {len(pypeit_file_models)} unique configs")
        for pypeit_file_model in pypeit_file_models:
            try:
                msgs.info(f"Creating tab {pypeit_file_model.config.name}")
                self.setUpdatesEnabled(False) # To prevent flickering when updating
                # Keep tabs sorted
                tab_index = None
                for i in range(len(self._file_tabs)):
                    if pypeit_file_model.config.name < self._file_tabs[i].name:
                        tab_index = i
                        break
                    elif pypeit_file_model.config.name == self._file_tabs[i].name:
                        msgs.bug(f"Duplicate name in tab list: {pypeit_file_model.config.name}")
                        return

                if tab_index == None:
                    tab_index = len(self._file_tabs)

                file_tab = PypeItFileTab(pypeit_file_model)
                pypeit_file_model.state_changed.connect(self.update_tab_status)

                self._file_tabs.insert(tab_index,file_tab)
                # Insert at tab_index + 1 because of the ObsLog tab
                self.insertTab(tab_index + 1, file_tab, file_tab.tab_name)

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
        for i in reversed(range(len(self._file_tabs))):
            if self._file_tabs[i].name in tab_list:
                # Remember to increment index by 1 because of the ObsLog tab
                self.removeTab(i+1)
                del self._file_tabs[i]
    
    def update_tab_status(self, name):
        """
        Update a tab's name to indicate if it has been changed but not saved.

        Args:
            tab_index (int): The index of the tab.
        """
        # Find the associated tab and set its name.
        for tab_index in range(len(self._file_tabs)):
            if self._file_tabs[tab_index].name == name:
                self.setTabText(tab_index+1, self._file_tabs[tab_index].tab_name)
                break

class SetupGUIMainWindow(QWidget):
    """Main window widget for the PypeIt Setup GUI

    Args:
        available_spectrographs (list of str): The list of spectrographs supported by PypeIt.
        controller (:class:`pypeit.setup_gui.controller.SetupGUIController`): The controller for the PypeITSetupGUI.
    """

    def __init__(self, model, controller):
        super().__init__(parent=None)

        self.layout = QVBoxLayout(self)
        self.model = model    
        self.controller = controller

        # Listen for operation_complete signals from the model, this
        # updates the progress dialog
        self.model.operation_complete.connect(self.operation_complete)

        # Layout the main window
        self.setup_view = PypeItSetupView(self.model) 
        self.layout.addWidget(self.setup_view)
        self.layout.addLayout(self._create_button_box())


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
        """Prompt the user if they want to save any unsaved data.

        Returns:
            DialogResponses: `SAVE` if the user wants to save, `ACCEPT` if they want to continue 
                             without saving, `CANCEL` if they want to cancel the current operation 
                             and to not lose any data.
        """
        response = QMessageBox.warning(self, self.tr("PypeIt Setup"), self.tr("There are unsaved changes to PypeIt setup files.\nDo you want to save then?"),
                                       QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel, QMessageBox.Save)
        if response == QMessageBox.Save:
            return DialogResponses.SAVE
        elif response == QMessageBox.Discard:
            return DialogResponses.ACCEPT
        else:
            return DialogResponses.CANCEL

    def prompt_for_open_file(self, caption, filter):
        """Opens a dialog to prompt the user for an existing file.
        
        Args:
            caption (str): A caption for the dialog
            filter (str):  A QFileDialog filter for a file to open. For example: "Text Files (`*`.txt)"

        Returns:
            str: The selected file, or None if the user pressed cancel.
        """
        # Get history from settings
        settings = QSettings()
        settings.beginGroup("OpenFile")
        history = settings.value("History")
        if history is None:
            history = []
        elif isinstance(history, str):
            history = [history]

        file_dialog = FileDialog(self, caption, QFileDialog.ExistingFile, filter=filter, history=history)        

        response =  file_dialog.show()
        if response != DialogResponses.CANCEL:
            # Add the selected dialog to history if it's new.
            # Note we only only wants paths in the history
            path = str(Path(file_dialog.selected_path).parent)
            if path not in history:
                history.append(path)
                settings.setValue("History", history)

        return file_dialog.selected_path


    def prompt_for_save_location(self, config_name, prompt_for_all=False):
        """Opens up a dialog requesting the user select a location
        to save a file. A history is maintained.

        Args:
            config_name (str): The name of the configuration being saved to a pypeit file.
            prompt_for_all (bool, Optional): Whether to prompt the user if this save location should
                                             apply to all subsequent files in the operation.

        Returns:
            tuple(str,DialogResponses): A tuple containing the selected location (None if the dialog was canceled).
                                        and the user's response: either CANCEL, ACCEPT, or ACCEPT_FOR_ALL.
        """

        # Get history
        settings = QSettings()
        settings.beginGroup("OutputDirectory")
        history = settings.value("History")
        if history is None:
            history = []
        elif not isinstance(history, list):
            history = [history]

        # Create the dialog.
        save_dialog = FileDialog(self, self.tr(f"Select location to save tab {config_name}."),
                                 QFileDialog.Directory, history=history, save_all=prompt_for_all)

        # Show the dialog.
        response = save_dialog.show()
        if response  != DialogResponses.CANCEL:            
            # the User selected a directory
            # save it to the history if it's new            
            if save_dialog.selected_path not in history:
                history.append(save_dialog.selected_path)
                settings.setValue("History", history)

        return (save_dialog.selected_path, response)

    def start_operation_progress(self, op_caption, max_progress_value, cancel_func):
        """Start displaying progress information for an operation. This uses the QProgressDialog, which will not
        display itself until a minimum amount of time has passed (currently 1s)
        
        Args:
            op_caption (str):         The name of the operation.
            max_progress_value (int): The maximum progress value (i.e. the value when done).
            cancel_func (:class:`collections.abc.Callable`):   A callable to deal with cancel being pressed in the 
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
        Increase the amount of progress being displayed in a progress dialog.

        Args:
            increase (int):          How much to increase the current progress by.
            message (str, Optional): A message indicating what step has been performed.
        """
        self.current_op_progress_dialog.setValue(self.current_op_progress_dialog.value() + increase)
        if message is not None:
            self.current_op_progress_dialog.setLabelText(message)

    def operation_complete(self):
        """
        Stop displaying progress for an operation because it has completed..
        """
        if self.current_op_progress_dialog is not None:
            self.current_op_progress_dialog.done(QDialog.Accepted)
            self.current_op_progress_dialog = None

    def update_save_tab_button(self):
        """Update the enabled/disabled state of the save tab button based on the
        current selected tab."""
        tab = self.setup_view.currentWidget()
        if tab.name == "ObsLog":
            self.saveTabButton.setEnabled(False)
        else:
            self.saveTabButton.setEnabled(tab.state == ModelState.CHANGED)

    def update_setup_button(self):
        """Enable/disable the setup button based on whether a spectrograph and raw
        data directories have been selected."""
        # Setup can only be run if the spectrograph is set and there's at least one
        # raw data directory
        if (self.model.spectrograph is not None and
            len(self.model.raw_data_directories) > 0):
            self.setupButton.setEnabled(True)
        else:
            self.setupButton.setEnabled(False)

    def update_buttons_from_model_state(self):
        """Update the enabled/disabled state of buttons based on the model state."""
        self.saveAllButton.setEnabled(self.setup_view.state==ModelState.CHANGED)
        self.clearButton.setEnabled(self.setup_view.state!=ModelState.NEW)
        self.update_save_tab_button()



    def _create_button_box(self):
        """Create the box with action buttons.
        
        Returns:
            QWidget: The widget with the action buttons for the GUI."""
            
        button_layout = QHBoxLayout()

        button = QPushButton(text = 'Open')
        button.setToolTip("Open a .pypeit file.")
        button.clicked.connect(self.controller.open_pypeit_file)
        button_layout.addWidget(button)
        self.openButton = button

        button = QPushButton(text = 'Clear')
        button.setToolTip("Clear everything and start with a blank slate.")
        button.setEnabled(False)
        button.clicked.connect(self.controller.clear)
        button_layout.addWidget(button)
        self.clearButton = button

        button = QPushButton(text = 'Run Setup')
        button.setToolTip("Scan the raw data and setup a .pypeit file for unique configuration.")
        button.setEnabled(False)
        button.clicked.connect(self.controller.run_setup)
        button_layout.addWidget(button)
        self.setupButton = button


        button = QPushButton(text = 'Save Tab')
        button.setToolTip("Save the curretly active tab.")
        button.setEnabled(False)
        button.clicked.connect(self.controller.save_one)
        button_layout.addWidget(button)
        self.saveTabButton = button

        button = QPushButton(text = 'Save All')
        button.setToolTip("Save all tabs with unsaved changes.")
        button.setEnabled(False)
        button.clicked.connect(self.controller.save_all)
        button_layout.addWidget(button)
        self.saveAllButton = button

        button_layout.addStretch()

        button = QPushButton(text = 'Exit')
        button.setToolTip("Quits this application.")
        button.clicked.connect(self.controller.exit)
        button_layout.addWidget(button)

        # Attach signals to enable/disable the Run Setup button and Save Tab button
        self.model.spectrograph_changed.connect(self.update_setup_button)
        self.model.raw_data_dirs_changed.connect(self.update_setup_button)
        self.setup_view.currentChanged.connect(self.update_save_tab_button)
        self.model.state_changed.connect(self.update_buttons_from_model_state)

        return button_layout

