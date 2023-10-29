import enum
from pathlib import Path
import traceback

from qtpy.QtWidgets import QGroupBox, QHBoxLayout, QVBoxLayout, QComboBox, QToolButton, QFileDialog, QWidget, QGridLayout, QFormLayout
from qtpy.QtWidgets import QMessageBox, QTabWidget, QTreeView, QLayout, QLabel, QScrollArea, QListView, QTableView, QPushButton, QProgressDialog, QDialog, QHeaderView, QSizePolicy, QCheckBox, QDialog
from qtpy.QtWidgets import QPlainTextEdit, QWidgetAction, QAction, QAbstractItemView, QStyledItemDelegate, QButtonGroup, QStyle, QTabBar
from qtpy.QtGui import QIcon, QKeySequence, QPalette, QColor, QValidator, QFont, QFontDatabase, QFontMetrics, QTextCharFormat, QTextCursor
from qtpy.QtCore import Qt, QSize, Signal,QSettings, QStringListModel, QAbstractItemModel, QModelIndex, QMargins, QSortFilterProxyModel, QRect

from pypeit.setup_gui.model import ModelState, PypeItMetadataModel, available_spectrographs
from pypeit import msgs

setup_gui_main_window = None

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
        save (bool, Optional):             Whether the dialog is a save dialog. If False it is treated as an open dialog. Defaults to False.
        ask_for_all (bool, Optional):      Whether the dialog should present a "Use this location for everything" option
                                           when saving. Defaults to False.
     """
    def __init__(self, parent, caption, file_mode, filter=None, history=[], save=False, ask_for_all=False):
        self._dialog = QFileDialog(parent, caption=parent.tr(caption), filter=filter)
        self._ask_for_all = ask_for_all

        self._dialog.setOption(QFileDialog.Option.DontUseNativeDialog)
        if save:
            self._dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
        else:
            self._dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptOpen)

        self._dialog.setFileMode(file_mode)
        if len(history) > 0:
            self._dialog.setDirectory(history[-1])
            self._dialog.setHistory(history)

        # Adding a checkbox like this only works for non-native dialogs. And if the 
        # default layout of the QFileDialog ever changes this could stop working.
        if self._ask_for_all:
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
            if self._ask_for_all and self.use_for_all_checkbox.isChecked():
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


class PathEditor(QWidget):
    """
    A custon widget that allows for entering a path. It has an editable combo box for entering 
    file locations, a browse button to use and a file dialog to enter the file location.
    The history of past locations is kept in the combo box list and file dialog history. 

    Args:
        parent(QWidget):                 The parent of this widget.
        browse_caption (str):            The caption text to use when searching for locations, also used as place holder
                                         text when no location has been set.
    """

    pathEntered = Signal(str)

    def __init__(self, browse_caption, parent=None):
        super().__init__(parent=parent)
        self.browse_caption = browse_caption

        # Setup the combo box
        layout = QHBoxLayout()
        layout.setSpacing(0)
        layout.setContentsMargins(0, 0, 0, 0)
        self._path = QComboBox(parent=self)
        self._path.setEditable(True)
        self._path.setInsertPolicy(QComboBox.InsertAtTop)

        # Setup history
        self._history = QStringListModel()
        self._path.setModel(self._history)
        if browse_caption is not None:
            self._path.lineEdit().setPlaceholderText(browse_caption)

        self._path.setCurrentIndex(-1)
        self._path.activated.connect(self._new_path)
        layout.addWidget(self._path)

        # Setup the browse button
        self._browse_button = QToolButton(parent=self)
        self._browse_button.setText(self.tr("Browse..."))
        self._browse_button.clicked.connect(self.browse)
        layout.addWidget(self._browse_button)

        self.setLayout(layout)
   
    def _add_path(self, new_path):
        """Set the path, adding it to history and signaling listeners.

        Args:
            new_path (str): The new path to add.
        """

        # Add to history if needed
        if new_path not in self._history.stringList():
            # It's not in the history, insert it into the 
            # combo box (which uses the history as its model)
            self._path.insertItem(0, new_path)

        self._path.setCurrentIndex(-1)
        self.pathEntered.emit(new_path)

    def setHistory(self, history):
        self._history = history
        self._path.setModel(history)
        self._path.setCurrentIndex(-1)

    def history(self):
        return self._history
    
    def browse(self):
        """Opens up a :class:`FileDialog` requesting the user select a location.

        Returns: 
            str: The new location, or None if Cancel was pressed.
        """
        browse_dialog = FileDialog(self, 
                                   caption = self.browse_caption, 
                                   file_mode=QFileDialog.Directory,
                                   history=self._history.stringList())

        if browse_dialog.show() == DialogResponses.ACCEPT:
            # the User selected a directory
            # Add it to our location list
            self._add_path(browse_dialog.selected_path)

            return browse_dialog.selected_path
        else:
            return None

    def _new_path(self, new_index):
        """
        Signal handler for when the combo box selects a new path.

        Args:
            new_index: The index within the combo box that was selected.
        """
        # Forward the signal to clients with the actual string value of the new
        # location
        if new_index != -1: 
            new_path = self._path.currentText()
            self._add_path(new_path)

class PypeItEnumListEditor(QWidget):
    def __init__(self, parent, allowed_values, num_lines=5):
        super().__init__(parent)
        self._values = set()
        self._allowed_values = allowed_values
        self._checkboxes = dict()
        layout = QVBoxLayout(self)
        layout.setSpacing(0)
        layout.setContentsMargins(0,0,0,0)
        # Wrap the list in a scroll area
        scroll_area = QScrollArea(parent=self)
        layout.addWidget(scroll_area)

        checkbox_container = QWidget()

        # Create the checkboxes for each allowable option
        self._button_group = QButtonGroup()
        self._button_group.setExclusive(False)
        checkbox_container.setAutoFillBackground(True)
        checkbox_layout=QVBoxLayout(checkbox_container)
        checkbox_layout.setContentsMargins(0,0,0,0)
        #self.setBackgroundRole(QPalette.ColorRole.Button)
        max_checkbox_width = 0
        for value in self._allowed_values:
            checkbox = QCheckBox(text=value, parent=checkbox_container)
            self._checkboxes[value] = checkbox
            self._button_group.addButton(checkbox)
            checkbox_layout.addWidget(checkbox)

            # Figure out the maximum size for a checkbox
            if checkbox.width() > max_checkbox_width:
                max_checkbox_width = checkbox.width()

        scroll_area.setWidget(checkbox_container)

        # Set the minimum height for this widget given requested # of lines
        # This assumes the checkboxes are the same height
        min_height = checkbox_layout.spacing()*(num_lines-1) + checkbox.height()*num_lines

        # Account for scrollbar            
        if scroll_area.horizontalScrollBar():
            min_height += scroll_area.horizontalScrollBar().sizeHint().height()

        # Account for margins
        min_height += scroll_area.contentsMargins().top() + scroll_area.contentsMargins().bottom()

        # Figure out the minimum width
        min_width = max_checkbox_width

        # Account for scrollbar size and margins
        if scroll_area.verticalScrollBar():
            min_width += scroll_area.verticalScrollBar().sizeHint().width()

        min_width += scroll_area.contentsMargins().left() + scroll_area.contentsMargins().right() 

        self.setFixedSize(min_width, min_height)
        msgs.info(f"Fixed size {min_width},{min_height}")
        self._button_group.buttonToggled.connect(self._choiceChecked)
        #msgs.info(f"My bg role/autofill: {self.backgroundRole()}/{self.foregroundRole()}/{self.autoFillBackground()} checkbox bg role/autofill: {checkbox.backgroundRole()}/{checkbox.foregroundRole()}/{checkbox.autoFillBackground()}")


    def _choiceChecked(self, widget, checked):
        value = widget.text()
        if checked:
            self._values.add(value)
        else:
            self._values.discard(value)


    def setSelectedValues(self, values):
        if values is None:
            self._values=set()
        else:
            self._values = set(values.split(","))
            for value in self._allowed_values:
                if value in self._values:
                    self._checkboxes[value].setChecked(True)
                else:
                    self._checkboxes[value].setChecked(False)

    def selectedValues(self):
        return ",".join(sorted(self._values))



class PypeItCustomEditorDelegate(QStyledItemDelegate):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def createEditor(self, parent,  option, index):
        
        model = index.model().sourceModel()
        
        column_name = model.getColumnNameFromNum(index)

        if column_name == "frametype":
            msgs.info("Creating enum list editor for frametype")
            return PypeItEnumListEditor(parent=parent, allowed_values=model.getAllFrameTypes())
        
        msgs.info(f"Creating default editor for {column_name}")
        return super().createEditor(parent, option, index)
    
    def setEditorData(self, editor, index):
        if isinstance(editor, PypeItEnumListEditor):
            msgs.info(f"Setting editor data {index.data(Qt.EditRole)}")
            editor.setSelectedValues(index.data(Qt.EditRole))
        else:
            msgs.info("Setting default editor data")
            super().setEditorData(editor, index)

    def setModelData(self,editor,model,index):
        if isinstance(editor, PypeItEnumListEditor):
            msgs.info(f"Setting choice model data: {editor.selectedValues()}")
            model.setData(index, editor.selectedValues())
        else:
            msgs.info("Setting default model data")
            super().setModelData(editor,model,index)

    def updateEditorGeometry(self, editor, option, index):
        if isinstance(editor, PypeItEnumListEditor):
            # The upper left coordinate of the editor depends on how well it fits
            # vertically in it's parent
            parent_geometry = editor.parent().geometry()
            editor_min_size = editor.minimumSize()

            msgs.info(f"Given rect: {(option.rect.x(), option.rect.y(), option.rect.width(), option.rect.height())}")
            msgs.info(f"parent_geometry: {(parent_geometry.x(), parent_geometry.y(), parent_geometry.width(), parent_geometry.height())}")
            msgs.info(f"editor min size: {editor_min_size.width()}, {editor_min_size.height()}")

            editor_x = option.rect.x()
            editor_y = option.rect.y()
            # Let the editor fill up the size of the cell if that's larger than
            # it's minimum width
            editor_width = max(editor_min_size.width(), option.rect.width())

            # Because the parent may be in a scroll area, the x,y could be negative,
            # set those to 0 so the editor doesn't appear outside the scroll area
            if editor_x < 0:
                editor_x = 0

            if editor_y < 0:
                editor_y = 0

            # Adjust the editor'x upper left corner so that the editor is vislbe for
            # cells along the bottom or right of the parent table
            if editor_x + editor_width > parent_geometry.bottomRight().x():
                editor_x = parent_geometry.bottomRight().x() - editor_width

            if editor_y + editor_min_size.height() > parent_geometry.bottomRight().y():
                editor_y = parent_geometry.bottomRight().y() - editor_min_size.height()

            geometry = QRect(editor_x, editor_y, editor_width, editor_min_size.height()) 
       
            msgs.info(f"Updating editor geometry to {(geometry.x(), geometry.y(), geometry.width(), geometry.height())}")
            editor.setGeometry(geometry)
        else:
            super().updateEditorGeometry(editor, option, index)


class PypeItMetadataView(QTableView):

    selectionUpdated = Signal()
    """Signal sent when items in the table are selected or deselected."""

    """QTableView to display file metadata.

    Args:
        parent (QtWidget):                                            The parent widget of the table view
        model  (:class:`pypeit.setup_gui.model.PypeItMetadataModel`): The model for the table. Optional, defaults to None.
    """
    def __init__(self, parent, model, controller):
        super().__init__(parent=parent)
        self._controller=controller
        self._controller.setView(self)
        self.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.setItemDelegate(PypeItCustomEditorDelegate(parent=self))
        self.setModel(model)
        self.addActions(controller.getActions(self))
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.ActionsContextMenu)

    def setModel(self, model):
        """Set the PypeItMetadataProxy model to use for the table.
        
        Args:
            model (:class:`pypeit.setup_gui.model.PypeItMetadataProxy`):  The model for the table.
        """
        # Use a sorting proxy model
        proxy_model = QSortFilterProxyModel()
        proxy_model.setSourceModel(model)
        sort_column = model.getColumnNumFromName('mjd')
        proxy_model.sort(sort_column, Qt.AscendingOrder)
        super().setModel(proxy_model)

        self.setSortingEnabled(True)
        self.horizontalHeader().setSortIndicator(sort_column, Qt.AscendingOrder)

    def selectionChanged(self, selected, deselected):
        """Event handler called by Qt when a selection change. Overriden from QTableView.
        
        Args:
            selected (QItemSelection): The items that are currently selected
            deselected (QItemSelection): The items that were deselected by the change.
        """
        super().selectionChanged(selected, deselected)
        self.selectionUpdated.emit()

    def selectedRows(self):
        """
        Return which rows in the PypeItMetadataPorxy model are selected.
        
        Return:
            list: list of rows indexes that are selected.
        """
        # Convert our selected indexes (whcih are in a proxy model) to the rows in our passed in model
        proxy_model = self.model()

        # Use a set to avoid adding duplicate rows caused by the selectedIndexes being per cell
        rows = set()
    
        for index in self.selectedIndexes():
            # The indexes are for each cell, but we're only interested in rows
            source_index = proxy_model.mapToSource(index)

            # Avoid adding duplicate rows caused by the selectedIndexes being per cell
            rows.add(source_index.row())
    
        return list(rows)

class ConfigValuesPanel(QGroupBox):

    """Scrollable panel to display configuration values for one or more frames.
    
    Args:
        spec_name (str):            Name of spectrograph for the configuration.
        config (list of tuple):     The name/value pairs for the configuration keys defined by the spectrograph.
        lines_to_display (int):     How many lines to display before scrolling.
        parent (QWidget, Optional): The parent widget, defaults to None.
    """
    def __init__(self, spec_name, config, lines_to_display, parent=None):

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
        form_widget_layout.addRow(self.tr("Spectrograph"), QLabel(spec_name))

        # We calculate our width to just fit around the text so we can stop the scroll area from
        # covering it with a scrollbar.
        fm = self.fontMetrics()
        max_key_width = fm.horizontalAdvance(self.tr("Spectrograph"))
        max_value_width = fm.horizontalAdvance(spec_name)

        # Add additional rows for configuration keys
        for key, value in config:
            label = QLabel(str(value))
            form_widget_layout.addRow(key, label)
            key_width = fm.horizontalAdvance(key)
            value_width = fm.horizontalAdvance(str(value))

            if key_width > max_key_width:
                max_key_width = key_width

            if value_width > max_value_width:
                max_value_width = value_width

        # Don't add extra margins in the FormLayout
        form_widget_layout.setContentsMargins(0, 0, 0, 0)

        scroll_area.setWidget(form_widget)

        # Set the minimum width  of the formwidget.        
        min_fw_width = form_widget_layout.horizontalSpacing()+max_key_width+max_value_width
        msgs.info(f"minWidth: {min_fw_width} max key width: {max_key_width} max_value width {max_value_width} horizontal spacing {form_widget_layout.horizontalSpacing()}")
    
        # Account for the scroll bar if needed
        if len(config) + 1 > lines_to_display:
            if scroll_area.verticalScrollBar():
                min_fw_width += scroll_area.verticalScrollBar().sizeHint().width()
                msgs.info(f"new minWidth: {min_fw_width} max key width: {max_key_width} max_value width {max_value_width} horizontal spacing {form_widget_layout.horizontalSpacing()}")

        form_widget.setMinimumWidth(min_fw_width)



        # Figure out the correct height for this panel, so that only the spectrograph and self.number_of_lines
        # config keys are visible

        # Find the minimum height of the form widget needed to hold the number of lines to display
        msgs.info(f"font height: {fm.height()} vertical spacing {form_widget_layout.verticalSpacing()}")
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

        # Set the minimum width of the form widget
        form_widget.setMinimumWidth(min_fw_width)

        # Set to fixed sizing policy
        policy = QSizePolicy()
        policy.setHorizontalPolicy(QSizePolicy.Minimum)
        policy.setVerticalPolicy(QSizePolicy.Fixed)
        policy.setControlType(QSizePolicy.DefaultType)
        self.setSizePolicy(policy)
        
        layout.addWidget(scroll_area)

class TabManagerBaseTab(QWidget):
    """Widget that acts as a tab for :class:`TabManagerWidget`. This defines the interface needed by TabManagerWidget
    and provides a default implementation of it.
    
    Args:

        state (:obj:`pypeit.setup_gui.model.ModelState`):  The state of the tab as displayed in the TabBar.
    """

    stateChanged = Signal(str, ModelState)
    """Signal(str): Signal sent when the state of the tab changes. The the name of the tab and its state is passed."""


    def __init__(self, parent=None, name="", closeable=False, state=ModelState.UNCHANGED):
        super().__init__(parent)
        self._name = name
        self._closeable=closeable
        self._state = state

    @property
    def name(self):
        return self._name

    @property
    def state(self):
        return self._state

    @property
    def closeable(self):
        return self._closeable

class PypeItFileView(TabManagerBaseTab):
    """Widget for displaying the information needed for one pypeit file. This includes
    the spectrograph, the configuration keys and values, the files that share that 
    configuration, and the PypeIt parameters.

    Args:
        (pypeit.setup_gui.model.PypeItFileModel): The model representing all the information needed for a .pypeit file.
        (pypeit.setup_gui.model.PypeItFileController): The controller for managing the user's interaction with a PypeIt file)
    """

    """Signal sent when items have been selected in the metadata table."""

    def __init__(self, model, controller):
        # Allow file tabs to be closed
        super().__init__(closeable=True)

        layout = QVBoxLayout(self) 
        self.model = model

        # Connect the model's state change signal to our own
        self.model.stateChanged.connect(self.stateChanged)

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

        # Add the ConfigValuesPanel, displaying the spectrograph + config keys.
        config_panel = ConfigValuesPanel(model.spec_name, model.config_values, 5, parent=self)
        third_row_layout.addWidget(config_panel)

        # Add the Raw Data directory panel to the third row, second column
        # This is not editable, because the user can add/remove directories by adding/removing individual
        # files in the metadata_file_table
        paths_group = QGroupBox(self.tr("Raw Data Directories"),self)
        paths_group_layout = QVBoxLayout(paths_group)        
        paths_viewer = QListView(paths_group)
        paths_viewer.setModel(model.paths_model)
        paths_viewer.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        paths_group_layout.addWidget(paths_viewer)   
        third_row_layout.addWidget(paths_group)

        # Make the paths wider than the config values panel
        third_row_layout.setStretch(1, 2)

        # Create a group box and table view for the file metadata table
        file_group = QGroupBox(self.tr("File Metadata"))
        file_group_layout = QVBoxLayout()
        self.file_metadata_table = PypeItMetadataView(self, model.metadata_model, controller.getMetadataController(model.metadata_model))
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

        self.model.stateChanged.connect(self.update_filename)


    """
    Signal handler that updates the filename when the underlying model is saved.
    """
    def update_filename(self):
        self.filename_value.setText(self.model.filename)

    @property
    def name(self):
        """
        str: The configuration name.
        """
        return self.model.name_stem

    @property
    def state(self):
        """
        :class:`pypeit.setup_gui.model.ModelState`): The state of this configuration's model. NEW, CHANGED, or UNCHANGED.
        """
        return self.model.state

    def selectedRows(self):
        """ Return the selected rows in the metadata table.

        Return:
            List of row indexes of the selected rows in the metadata table.
        """
        return self.file_metadata_table.selectedRows()

class ObsLogView(TabManagerBaseTab):
    """Widget for displaying the observation log for raw data files for the same spectrograph but 
    potentially different observing configurations.

    Args:
        model (:class:`pypeit.setup_gui.model.PypeItObsLogModel`): Model object for a PypeIt Setup GUI.
        controller (:class:`pypeit.setup_gui.model.PypeItObsLogController`): Controller object for the PypeIt Setup GUI.
        parent (QWidget): The parent widget of the tab.
    """

    """Signal sent when items are selected in the metadata table."""

    def __init__(self, model, controller, parent=None):
        super().__init__(parent=parent,name="ObsLog")

        self._controller = controller

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

        # Create a Group Box to group the paths editor and viewer
        paths_group = QGroupBox(self.tr("Raw Data Directories"),self)
        paths_group_layout = QVBoxLayout(paths_group)        

        self.paths_editor = PathEditor(browse_caption=self.tr("Choose raw data directory"), parent=self)
        self.paths_editor.setHistory(PersistentStringListModel("RawDataDirectory", "History"))
        self.paths_editor.setEnabled(False)
        self.paths_editor.pathEntered.connect(controller.addNewPath)

        paths_group_layout.addWidget(self.paths_editor)


        self._paths_viewer = QListView(paths_group)
        self._paths_viewer.setModel(model.paths_model)
        fm = self.fontMetrics()
        # Only display 5 paths
        lines =5
        self._paths_viewer.setFixedHeight(fm.height()*lines+self._paths_viewer.spacing()*(lines-1))
        paths_group_layout.addWidget(self._paths_viewer)

        # Add action for removing a path
        self._paths_viewer.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        delete_action = QAction(self._paths_viewer, text=self.tr("Delete selected"))
        delete_action.setToolTip(self.tr("Delete selected path"))
        delete_action.setShortcut(QKeySequence.StandardKey.Delete)
        delete_action.triggered.connect(self._deletePaths)
        self._paths_viewer.addAction(delete_action)
        self._paths_viewer.setContextMenuPolicy(Qt.ContextMenuPolicy.ActionsContextMenu)


        # Add the Raw Data directory panel to the first row, second column
        top_row_layout.addWidget(paths_group)

        # Make the metadata wider than the spectrograph
        top_row_layout.setStretch(1, 2)


        # Add the File Metadata box in the second row
        file_group = QGroupBox(self.tr("File Metadata"))
        file_group_layout = QHBoxLayout()        
        self.obslog_table = PypeItMetadataView(file_group, model.metadata_model, controller.getMetadataController(model.metadata_model))
        file_group_layout.addWidget(self.obslog_table)
        file_group.setLayout(file_group_layout)
        layout.addWidget(file_group)
        # Make File Metadata taller than the spectrograph/raw data paths row
        layout.setStretch(1,4)
        self.setModel(model)

        # Update model with new spectrograph/data paths
        self.spectrograph.textActivated.connect(controller.setSpectrograph)
        self.spectrograph.textActivated.connect(self.update_raw_data_paths_state)

    def _deletePaths(self, parent):
        msgs.info(f"Delete selection")
        selection = self._paths_viewer.selectedIndexes()
        rows = [index.row() for index in selection]
        self._controller.removePaths(rows)

    def setModel(self,model):
        """Set a new model for file metadata.

        Args:
            model (:class:`pypeit.setup_gui.model.PypeItSetupModel`): The new metadata model
        """
        self.model=model
        if model.spec_name is not None:
            self.spectrograph.setCurrentIndex(self.spectrograph.findText(model.spec_name))
            msgs.info(f"Set current text to {model.spec_name}, current index {self.spectrograph.currentIndex()}")
            self.update_raw_data_paths_state()
        self.obslog_table.setModel(model.metadata_model)
        self._controller.setModel(model)
        # Update based on changes to the model
        model.spectrograph_changed.connect(self.update_from_model)
        model.paths_model.dataChanged.connect(self.update_from_model)
        model.paths_model.rowsRemoved.connect(self.update_from_model)
        model.paths_model.modelReset.connect(self.update_from_model)


    def update_raw_data_paths_state(self):
        """Enable/Disable the raw data paths location panel based on the model state. """
        if self.model.state == ModelState.NEW:
            if self.spectrograph.currentIndex() == -1:
                self.paths_editor.setEnabled(False)
            else:
                self.paths_editor.setEnabled(True)
        else:
            self.paths_editor.setEnabled(True)

    def update_from_model(self):
        """
        Updates the spectrograph and raw data location widgets based on model updates.
        """
        # Update the current spectrograph
        if self.model.spec_name is None:
            self.spectrograph.setCurrentIndex(-1)
        else:
            self.spectrograph.setCurrentText(self.model.spec_name)

        # Disable changing the spectrograph if the model isn't in the new state
        if self.model.state == ModelState.NEW:
            self.spectrograph.setEnabled(True)
        else:
            self.spectrograph.setEnabled(False)

        # Get the raw data paths state
        self.update_raw_data_paths_state()

    def selectedRows(self):
        """ Return the selected rows in the metadata table.

        Return:
            List of row indexes of the selected rows in the metadata table.
        """
        # Return the metadata view's selected items
        return self.obslog_table.selectedRows()

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
                return QValidator.Acceptable, str_input, int_input
            else:
                for spectrograph in self._supported_spectrographs:
                    if spectrograph.startswith(str_input.lower()):
                        return QValidator.Intermediate, str_input, int_input
        return QValidator.Invalid, str_input, int_input

class TabManagerWidget(QTabWidget):
    """
    Widget which manages tabs for the obslog and pypeit files.
    It extends the QTabWidget functionality by allowing tabs by
    adding widgets to add and remove tabs as well as specifying
    the modification state with a "*" when a tab hasn't been saved.

    Args:
        model (:class:`pypeit.setup_gui.model.SetupGUIStateModel`): Model object for a PypeIt Setup GUI.
        controller (:class:`pypeit.setup_gui.model.SetupGUIControllerl`): Controller object for the PypeIt Setup GUI.
    """

    tabCreateRequest = Signal()

    def __init__(self, parent, tab_position=QTabWidget.South):
        super().__init__(parent=parent)
        
        self._tabNames = []
        self.setTabPosition(tab_position)

        # Allow tabs to be closed
        self.setTabsClosable(True)

        # Add a + tab to create new tabs
        self.addNewTab(TabManagerBaseTab(parent=self,name="+"))
        self.tabBarClicked.connect(self.checkNewTabClicked)
        self.currentChanged.connect(self.checkIfNewTabCurrent)

    def checkNewTabClicked(self, index):
        if index == self.count() - 1:
            # Create as new tab. The new tab model will send the signals needed
            # to create the view and add it
            self.tabCreateRequest.emit()

    def checkIfNewTabCurrent(self, index):
        # Try to prevent the + tab from being visible
        if self.count() > 1:
            if self.widget(index).name == "+":
                # Make the first tab visible instead
                self.setCurrentIndex(0)

    def addNewTab(self, tab):
        # Insert before the "+" tab. If this is the first tab being entered
        # (ie the + tab), it is appended
        index = self.count()-1
        index=self.insertTab(index, tab, tab.name)
        msgs.info(f"Added {tab.name} at index {index}")
        self._tabNames.insert(index,tab.name)
        self.updateTabText(tab.name,tab.state)
        if tab.closeable:
            tab.stateChanged.connect(self.updateTabText)
        else:
            self._hideTabCloseButton(index, always=True)
        return index

    def closeTab(self, tab_name):
        index = self._tabNames.index(tab_name)
        tab = self.widget(index)
        if tab.closeable:
            tab.stateChanged.disconnect(self.updateTabText)
        self.removeTab(index)
        del self._tabNames[index]

    def updateTabText(self, tab_name, tab_state):
        index = self._tabNames.index(tab_name)
        tab = self.widget(index)
        if tab_state != ModelState.UNCHANGED:
            self.setTabText(index, "*" + tab.name)
        else:
            self.setTabText(index, tab.name)

    def setNewTabsEnabled(self, enabled):
        """Sets whether or not new tabs can be added to the tab widget.
        This will enable/disable the "+" tab.
        Args:
            enabled (bool): True if new tabs can be added, Fale otherwise.
        """
        if self.count() > 0:
            self.setTabEnabled(self.count()-1,enabled)

    def _hideTabCloseButton(self, tab_pos, always=False):
        """Hides the tab close button on a tab. Used for the ObsLog tab.
        
        Args:
            tab_pos (int): The index of the tab.
            always (bool): True if the close tab button should always be hidden. If True
                           the tab will not reserve any space for the button to appear.
        """
        close_button_positon = QTabBar.ButtonPosition(self.style().styleHint(QStyle.SH_TabBar_CloseButtonPosition))
        close_button = self.tabBar().tabButton(tab_pos, close_button_positon)
        if always:
            close_button.setFixedWidth(tab_pos)
        close_button.hide()



class LogWindow(QWidget):
    """Window showing PypeIt log messages as they occur. The log messages can also be saved to a text file.
    
    Args:
        main_window(SetupGUIMainWindow): The main window of the PypeIt Setup GUI.
        logBuffer(:obj:`pypeit.setup_gui.model.LogBuffer`): The log buffer holding the log messages.
    """

    closed = Signal()
    """Signal sent when the log window is closed."""

    newMessage = Signal(str)
    """Signal sent when a new log message is logged."""

    def __init__(self, main_window, logBuffer):
        super().__init__()
        self._logBuffer = logBuffer
        self._main_window = main_window
        self.setWindowIcon(QIcon(str(Path(__file__).parent / "images/window_icon.png")))
        self.setWindowTitle(self.tr("PypeIt Setup Log"))
        logLayout = QVBoxLayout(self)

        # Get a fixed width font
        fixed_font = QFont()
        fixed_font.setFamilies(["Monospace", "Courier"])
        fixed_font.setFixedPitch(True)
        fixed_font.setPointSize(12)

        # Create the log viewer widget
        self.logViewer=QPlainTextEdit(parent=self)
        self.logViewer.setFont(fixed_font)
        self.logViewer.document().setMaximumBlockCount(self._logBuffer.maxlen)
        self.logViewer.setReadOnly(True)

        # Set the log window size to 100x50 characters
        viewer_margins = self.logViewer.contentsMargins()
        layout_margins = logLayout.contentsMargins()
        parent_margins = self.contentsMargins()
        font_metrics = self.logViewer.fontMetrics()
        char_width = font_metrics.averageCharWidth()*100
        char_height = font_metrics.height()*50
        
        self.resize(parent_margins.left() + parent_margins.right() +
                    layout_margins.left() + layout_margins.right() +
                    viewer_margins.left() + viewer_margins.right() +
                    char_width, 
                    parent_margins.top() + parent_margins.bottom() +
                    layout_margins.top() + layout_margins.bottom() +
                    viewer_margins.top() + viewer_margins.bottom() +
                    char_height)

        logLayout.addWidget(self.logViewer)

        # Listen for new log messages and append them to the viewer
        # This is done with a callbackto self._messageLogged, which then
        # emits an event from self.newMessage to notify self._addMessage.
        # This is done so that log messages from a different thread are
        # properly queued to the GUIs event thread.
        self.newMessage.connect(self._addMessage,Qt.QueuedConnection)
        self._logBuffer.watch("log window", None, self._messageLogged)

        # Fill the log viewer with previously logged messages.
        for message in self._logBuffer:
            self._addMessage(message)


        # Buttons
        buttonLayout = QHBoxLayout()
        logLayout.addLayout(buttonLayout)

        # Button to save the log
        self.saveButton = QPushButton(text=self.tr("Save Log..."))
        self.saveButton.clicked.connect(self.save)
        buttonLayout.addWidget(self.saveButton)
        buttonLayout.addStretch()

        # Close button
        self.closeButton = QPushButton(text=self.tr("Close"))
        self.closeButton.clicked.connect(self.close)
        buttonLayout.addWidget(self.closeButton)

    def save(self):
        """Save the log to a file."""

        # Get history
        settings = QSettings()
        settings.beginGroup("LogDirectory")
        history = settings.value("History")
        if history is None:
            history = []
        elif not isinstance(history, list):
            history = [history]

        # Create the dialog.
        save_dialog = FileDialog(self, self.tr(f"Enter name for logfile."),
                                 QFileDialog.AnyFile, history=history, filter="Log Files (`*`.log)", save=True)

        # Show the dialog.
        response = save_dialog.show()
        if response  != DialogResponses.CANCEL:            
            # the User selected a directory
            # save it to the history if it's new            
            if save_dialog.selected_path not in history:
                history.append(save_dialog.selected_path)
                settings.setValue("History", history)

            # Save the file
            try:
                with open(save_dialog.selected_path, "w") as f:
                    for message in self._logBuffer:
                        f.write(message)
                msgs.info(f"Log saved to {save_dialog.selected_path}.")
            except Exception as e:
                self._main_window.display_error(str(e))

    def closeEvent(self, event):
        """Event handle for closing the log window. Overridden from QWidget.

        Args:
            event (QEvent): The close event. Not used by this implementation.
        """
        # Stop watching the log, and notify the main window.
        self.newMessage.disconnect(self._addMessage)
        self._logBuffer.unwatch("log window")
        self.closed.emit()
        return super().closeEvent(event)

    def _messageLogged(self, message):
        """Callback that is notified by the log buffer for each message logged.

        Args:
            message (str): The log message
        """
        # Emit the message so it can be picked up in the GUI's event thread.
        self.newMessage.emit(message)

    def _addMessage(self, message):
        """Signal handler that is notified on the GUI event thread of a new message.
        
        Args:
            message (str): The message logged.
        """

        # Add the message to our internal text document
        self.logViewer.appendPlainText(message.strip())


class SetupGUIMainWindow(QWidget):
    """Main window widget for the PypeIt Setup GUI

    Args:
        model (:class:`pypeit.setup_gui.model.SetupGUIStateModel`): The model for the PypeitSetupGUI.
        controller (:class:`pypeit.setup_gui.controller.SetupGUIController`): The controller for the PypeitSetupGUI.
    """

    def __init__(self, model, controller):
        super().__init__(parent=None)

        self.layout = QVBoxLayout(self)
        self.model = model    
        self.controller = controller

        # Create the initial observation log tab
        self._obs_log_tab = ObsLogView(model = model.obslog_model,
                                       controller = controller.getObsLogController(model.obslog_model))

        # Layout the main window
        self.tab_widget = TabManagerWidget(parent=self)
        index = self.tab_widget.addNewTab(self._obs_log_tab)
        self.tab_widget.setCurrentIndex(index)
        self.tab_widget.tabCreateRequest.connect(self.controller.createNewPypeItFile)

        # enable/disable adding new files based on whether the spectrograph has been set
        self.model.obslog_model.spectrograph_changed.connect(self.update_new_file_allowed)
        self.update_new_file_allowed()

        self.layout.addWidget(self.tab_widget)

        # Monitor the current tab
        self._current_tab = self._obs_log_tab

        # Create the row of buttons for user actions
        self.layout.addLayout(self._create_button_box())

        # For viewing the log
        self._logWindow = None

        # For displaying operation progress
        self.current_op_progress_dialog = None

        # Monitor the model for new or closed files, and update the tab widget accordingly
        # These are queued because this can happen in a background operation thread
        self.model.filesAdded.connect(self.create_file_tabs)
        self.model.filesDeleted.connect(self.delete_tabs)

        # Setup application/window icon TODO this doesn't work in windows. Mac???
        self.setWindowIcon(QIcon(str(Path(__file__).parent / "images/window_icon.png")))
        self.setWindowTitle(self.tr("PypeIt Setup"))

        self.resize(1650,900)

    def update_new_file_allowed(self):
        """Signal handler to enable/disable adding a new file based on whether the spectrograph has been set."""        
        self.tab_widget.setNewTabsEnabled(self.model.obslog_model.spec_name is not None)
        

    def create_progress_dialog(self, op_caption, max_progress_value, cancel_func):
        """Start displaying progress information for an operation. This uses the QProgressDialog, which will not
        display itself until a minimum amount of time has passed (currently 1s)
        
        Args:
            op_caption (str):         The name of the operation.
            max_progress_value (int): The maximum progress value (i.e. the value when done).
            cancel_func (:class:`collections.abc.Callable`):   A callable to deal with cancel being pressed in the 
                                                               progress dialog.
        """
        msgs.info(f"Starting operation {op_caption} max progress: {max_progress_value}")
        self.current_op_progress_dialog = QProgressDialog(self.tr(op_caption), self.tr("Cancel"), 0, max_progress_value, parent=self)
        self.current_op_progress_dialog.setMinimumWidth(380)
        self.current_op_progress_dialog.setWindowTitle(op_caption)
        self.current_op_progress_dialog.setMinimumDuration(1000)
        self.current_op_progress_dialog.setValue(0)            
        self.current_op_progress_dialog.canceled.connect(cancel_func)

    def show_operation_progress(self, increase, message=None):
        """
        Increase the amount of progress being displayed in a progress dialog.

        Args:
            increase (int):          How much to increase the current progress by.
            message (str, Optional): A message indicating what step has been performed.
        """
        msgs.info(f"dialog is none {self.current_op_progress_dialog is None}")
        if self.current_op_progress_dialog is not None:
            msgs.info(f"increase {increase} message{message} current value {self.current_op_progress_dialog.value()}")
            self.current_op_progress_dialog.setValue(self.current_op_progress_dialog.value() + increase)
            if message is not None:
                self.current_op_progress_dialog.setLabelText(message)

    def operation_complete(self):
        """
        Stop displaying progress for an operation because it has completed..
        """
        msgs.info(f"Ending operation, dialog is none {self.current_op_progress_dialog is None}")
        if self.current_op_progress_dialog is not None:
            self.current_op_progress_dialog.done(QDialog.Accepted)
            self.current_op_progress_dialog = None

    def update_save_tab_button(self):
        """Update the enabled/disabled state of the save tab button based on the
        current selected tab."""
        tab = self.tab_widget.currentWidget()
        if tab.name == "ObsLog":
            self.saveTabButton.setEnabled(False)
        else:
            self.saveTabButton.setEnabled(tab.state != ModelState.UNCHANGED)

    def update_setup_button(self):
        """Enable/disable the setup button based on whether a spectrograph and raw
        data directories have been selected."""
        # Setup can only be run if the spectrograph is set and there's at least one
        # raw data directory
        msgs.info(f"Checking setup button status spec: {self.model.obslog_model.spec_name} dirs {self.model.obslog_model.raw_data_directories}")
        if (self.model.obslog_model.spec_name is not None and
            len(self.model.obslog_model.raw_data_directories) > 0):
            self.setupButton.setEnabled(True)
        else:
            self.setupButton.setEnabled(False)

    def update_buttons_from_model_state(self):
        """Update the enabled/disabled state of buttons based on the model state."""
        self.saveAllButton.setEnabled(self.model.state==ModelState.CHANGED)
        self.clearButton.setEnabled(self.model.state!=ModelState.NEW)
        self.update_save_tab_button()

    def _showLog(self):
        """Signal handler that opens the log window."""
        if self._logWindow is not None:
            self._logWindow.activateWindow()
            self._logWindow.raise_()
        else:
            self._logWindow = LogWindow(self, self.model.log_buffer)
            self._logWindow.closed.connect(self._logClosed)
            self._logWindow.show()

    def _logClosed(self):
        """Signal handler that clears the log window when it closes."""
        self._logWindow = None
            
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

        button = QPushButton(text = 'View log')
        button.setToolTip("Opens a window containing the log.")
        button.clicked.connect(self._showLog)
        button_layout.addWidget(button)
        self.logButton = button

        button = QPushButton(text = 'Exit')
        button.setToolTip("Quits this application.")
        button.clicked.connect(self.controller.exit)
        button_layout.addWidget(button)

        # Monitor when new files are added and removed,
        # so we can update buttons.
        self.model.obslog_model.paths_model.dataChanged.connect(self.update_setup_button)
        self.model.obslog_model.paths_model.rowsRemoved.connect(self.update_setup_button)
        self.model.obslog_model.paths_model.modelReset.connect(self.update_setup_button)

        # Monitor spectrograph changes
        self.model.obslog_model.spectrograph_changed.connect(self.update_setup_button)

        # Monitor the current tab to update the save_tab button
        self.tab_widget.currentChanged.connect(self.update_save_tab_button)

        # Monitor the application's NEW/CHANGED/UNCHANGED state to enable/disable buttons
        self.model.stateChanged.connect(self.update_buttons_from_model_state)

        # Setup may have already been run from the command line, so update button status on init
        self.update_setup_button()
        self.update_buttons_from_model_state()
        return button_layout

    def create_file_tabs(self, pypeit_file_models):
        """
        Create new tabs for new unique configurations found either by openeing a pypeit file or
        by reading raw data directories.

        Args:
            pypeit_file_models (list of :class:`pypeit.setup_gui.model.PypeItFileModel`): Models for the tabs to add.
        """
        msgs.info(f"create_file_tabs for {len(pypeit_file_models)} unique configs")
        try:
            self.tab_widget.setUpdatesEnabled(False) # To prevent flickering when updating
            for model in pypeit_file_models:
                new_tab =  PypeItFileView(model, self.controller.getPypeItFileController(model))
                self.tab_widget.addNewTab(new_tab)
        finally:
            self.tab_widget.setUpdatesEnabled(True) # To allow redrawing after updating

    def delete_tabs(self, file_list):
        """
        Delete any tabs no longer in the model.

        Args:
            tab_list (list of str): List of the configuration names removed.
        """
        msgs.info(f"View Deleting tabs {file_list}")
        if len(file_list) == 0:
            return

        for file_name in  file_list:
            self.tab_widget.closeTab(file_name)

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
        else:
            file_dialog.selected_path = None

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
                                 QFileDialog.Directory, history=history, save=True, ask_for_all=prompt_for_all)

        # Show the dialog.
        response = save_dialog.show()
        if response  != DialogResponses.CANCEL:            
            # the User selected a directory
            # save it to the history if it's new            
            if save_dialog.selected_path not in history:
                history.append(save_dialog.selected_path)
                settings.setValue("History", history)
        else:
            save_dialog.selected_path = None

        return (save_dialog.selected_path, response)

