"""
The view portion of the PypeIt Setup GUI.  Responsible for displaying information to the user, and forwarding user input to the controller.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pathlib import Path

from qtpy.QtWidgets import QGroupBox, QHBoxLayout, QVBoxLayout, QComboBox, QToolButton, QFileDialog, QWidget, QGridLayout, QFormLayout
from qtpy.QtWidgets import QMessageBox, QTabWidget, QTreeView, QLayout, QLabel, QScrollArea, QListView, QTableView, QPushButton, QStyleOptionButton, QProgressDialog, QDialog, QHeaderView, QSizePolicy, QCheckBox, QDialog
from qtpy.QtWidgets import QAction, QAbstractItemView, QStyledItemDelegate, QButtonGroup, QStyle, QTabBar,QAbstractItemDelegate
from qtpy.QtGui import QIcon,QMouseEvent, QKeySequence, QPalette, QColor, QValidator, QFont, QFontDatabase, QFontMetrics, QTextCharFormat, QTextCursor
from qtpy.QtCore import Qt, QObject, QSize, Signal,QSettings, QStringListModel, QAbstractItemModel, QModelIndex, QMargins, QSortFilterProxyModel, QRect

from pypeit.spectrographs import  available_spectrographs

from pypeit.setup_gui.model import ModelState, PypeItMetadataModel
from pypeit.setup_gui.text_viewer import LogWindow, TextViewerWindow
from pypeit.setup_gui.dialog_helpers import DialogResponses, FileDialog, PersistentStringListModel
from pypeit import msgs

def debugSizeStuff(widget:QWidget, name="widget"):
    """Helper method for logging sizxing information about a wdiget and its layout."""
    msgs.info(f"{name} (width/height): {widget.width()}/{widget.height()} geometry x/y/w/h: {widget.geometry().x()}/{widget.geometry().y()}/{widget.geometry().width()}/{widget.geometry().height()} min w/h {widget.minimumWidth()}/{widget.minimumHeight()} hint w/h {widget.sizeHint().width()}/{widget.sizeHint().height()} min hint w/h {widget.minimumSizeHint().width()}/{widget.minimumSizeHint().height()} cm tlbr: {widget.contentsMargins().top()}/{widget.contentsMargins().left()}/{widget.contentsMargins().bottom()}/{widget.contentsMargins().right()} frame w/h {widget.frameSize().width()}/{widget.frameSize().height()}")
    layout = widget.layout()
    if layout is None:
        msgs.info(f"{name} layout is None")
    else:
        msgs.info(f"{name} layout size constraint {layout.sizeConstraint()} spacing: {layout.spacing()} cm: tlbr {layout.contentsMargins().top()}/{layout.contentsMargins().left()}/{layout.contentsMargins().bottom()}/{layout.contentsMargins().right()} totalMinSize (w/h): {layout.totalMinimumSize().width()}/{layout.totalMinimumSize().width()} totalMaxSize (w/h): {layout.totalMaximumSize().width()}/{layout.totalMaximumSize().width()} totalHint (w/h): {layout.totalSizeHint().width()}/{layout.totalSizeHint().width()}")

def calculateButtonMinSize(button_widget : QPushButton) -> QSize:
    """Calculates and sets the minimum size of a budget widget
    
    Qt has code in QCommonStyle to set this size for a button, but I kept discovering that it would report a much
    larger size for some reason. So this method exists to fix that.
    
    Args:
        button_widget: The button to set the minimum size for. It should already have it's text set.
    
    Return:
        The minimum size that was calcualted for the button
    """
    # Get the size of the button's text
    fm = button_widget.fontMetrics()
    text_size = fm.size(Qt.TextFlag.TextShowMnemonic,button_widget.text())

    # Get the style sizes for the frame, margin, and (if applicable) the default indicator. These are integer values
    style_options = QStyleOptionButton()
    style_options = button_widget.initStyleOption(style_options)
    button_margin = button_widget.style().pixelMetric(QStyle.PixelMetric.PM_ButtonMargin, style_options,button_widget)
    button_default_frame = button_widget.style().pixelMetric(QStyle.PixelMetric.PM_DefaultFrameWidth,style_options,button_widget)
    if button_widget.isDefault():
        default_indicator = button_widget.style().pixelMetric(QStyle.PixelMetric.PM_ButtonDefaultIndicator,style_options,button_widget)
    else:
        default_indicator = 0


    # The QT code doubles the frame size but not the margin, so we do the same
    min_size = QSize(text_size.width() + button_margin + button_default_frame*2 + default_indicator*2,
                     text_size.height() + button_margin + button_default_frame*2 + default_indicator*2)
    msgs.info(f"Calculated button {button_widget.text()} minimum size ({min_size.width()}/{min_size.height()}) with text_size ({text_size.width()}/{text_size.height()}) margin size ({button_margin}) frame width ({button_default_frame}) and default indicator width ({default_indicator})")
    
    return min_size
    


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
    """Signal sent when a path has been added."""

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
        """Sets the past history of the PathEditor widget.
        
        Args:
            history (QStringListModel): A string list model containing the history of the widget."""
        self._history = history
        self._path.setModel(history)
        self._path.setCurrentIndex(-1)

    def history(self):
        """Returns the past history of the PathEditor widget.
        
        Returns:
            QStringListModel: A stringh list model with the history of the widget.
        """
        return self._history
    
    def browse(self):
        """Opens up a :class:`FileDialog` requesting the user select a location.

        Returns: 
            str: The new location, or None if Cancel was pressed.
        """
        browse_dialog = FileDialog(self, 
                                   caption = self.browse_caption, 
                                   file_mode=QFileDialog.Directory,
                                   history=self._history)

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
    """Widget for editing a enumerated list of values by checking the values
    on or off with a checkbox.
    
    Args:
        parent (QWidget): 
            The parent of the editor.
        allowed_values (list of str):
            The list of allowed values in the enumeration.
        index (QModelIndex):
            The index of the item being edited.
        num_lines (int):
            The number of lines to display. Any other lines
            will be reachable by scrolling.
    """

    closed = Signal(QWidget, bool)
    """
    Signal sent when the user closes the editor with the OK or CANCEL button.
    The signal will provide the editor that was closed and a boolean that will
    be True if the change was accepted or False if it was canceled.
    """

    def __init__(self, parent, allowed_values, index, num_lines):
        super().__init__(parent)
        self.index=index
        self._values = set()
        self._allowed_values = allowed_values
        self._checkboxes = dict()
        self.setBackgroundRole(QPalette.ColorRole.Window)
        self.setAutoFillBackground(True)
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

        max_checkbox_width = 0
        for value in self._allowed_values:
            checkbox = QCheckBox(text=value, parent=checkbox_container)
            self._checkboxes[value] = checkbox
            self._button_group.addButton(checkbox)
            checkbox_layout.addWidget(checkbox)

            # Figure out the maximum size for a checkbox
            if checkbox.width() > max_checkbox_width:
                max_checkbox_width = checkbox.width()

        msgs.info(f"Max checkbox width: {max_checkbox_width}")
        scroll_area.setWidget(checkbox_container)

        # Figure out the minimum width
        min_width = max_checkbox_width

        # Account for scrollbar size and margins
        if scroll_area.verticalScrollBar():
            min_width += scroll_area.verticalScrollBar().sizeHint().width()

        min_width += scroll_area.contentsMargins().left() + scroll_area.contentsMargins().right() 

        # Add Ok and cancel buttons at the bottom, outside the scroll area
        ok_cancel_container = QWidget(parent=checkbox_container)
        ok_cancel_layout = QHBoxLayout()

        accept_button=QPushButton(text="OK")
        accept_button.setDefault(True)
        accept_button.clicked.connect(self._accepted)
        cancel_button=QPushButton(text="Cancel")
        cancel_button.clicked.connect(self._canceled)
        ok_cancel_layout.addWidget(accept_button)
        ok_cancel_layout.addWidget(cancel_button)
        ok_cancel_container.setLayout(ok_cancel_layout)

        ok_cancel_layout_margins = ok_cancel_layout.contentsMargins()

        # Use small margins along the left/right
        ok_cancel_layout_margins.setRight(1)
        ok_cancel_layout_margins.setLeft(1)
        ok_cancel_layout.setContentsMargins(ok_cancel_layout_margins)

        # Make sure the minimum width doesn't truncate the button's text
        ok_button_min_size = calculateButtonMinSize(accept_button)
        cancel_button_min_size = calculateButtonMinSize(cancel_button)
        button_min_width = max(ok_button_min_size.width(), cancel_button_min_size.width())    
    
        ok_cancel_container_min_width = button_min_width*2 + ok_cancel_layout.spacing() + ok_cancel_layout_margins.left() + ok_cancel_layout_margins.right()
        msgs.info(f"Okay cancel container min_width: {ok_cancel_container_min_width}")
        if min_width < ok_cancel_container_min_width:
            min_width = ok_cancel_container_min_width

        # Set the okay/cancel button container's min width to keep Qt
        # from making the buttons bigger than they need to be.
        ok_cancel_container.setMinimumWidth(min_width)
        layout.addWidget(ok_cancel_container)

        # Set the minimum height for this widget given requested # of lines
        # This assumes the checkboxes are the same height
        min_height = checkbox_layout.spacing()*(num_lines-1) + checkbox.height()*num_lines

        # Account for margins
        min_height += scroll_area.contentsMargins().top() + scroll_area.contentsMargins().bottom()

        # Account for buttons
        min_height += ok_button_min_size.height() + ok_cancel_layout.contentsMargins().top() + ok_cancel_layout.contentsMargins().bottom()

        self.setMinimumSize(min_width, min_height)
        self._button_group.buttonToggled.connect(self._choiceChecked)

        msgs.info(f"min_width/height: {min_width}/{min_height}")
        debugSizeStuff(self, "Enum Editor")
        debugSizeStuff(checkbox_container, "Checkbox Container")
        debugSizeStuff(ok_cancel_container, "OK/Cancel Container")
        debugSizeStuff(accept_button, "OK Button")
        debugSizeStuff(cancel_button, "Cancel Button")
        

    def _accepted(self, *args):
        """Signal handler for when the "OK" button is clicked."""
        self.closed.emit(self, True)

    def _canceled(self, *args):
        """Signal handler for when the "Cancel" button is clicked."""
        self.closed.emit(self, False)

    def _choiceChecked(self, widget, checked):
        """Signal handler for when one of the enumerated values is checked on or off.
        Args:
            widget (QCheckBox):
                The widget that checked or unchecked.
            checked (bool):
                 True if the widget is checked, false if it is not.
        """
        value = widget.text()
        if checked:
            self._values.add(value)
        else:
            self._values.discard(value)

    def setSelectedValues(self, values):
        """Set what values of the enumeration are selected.
        
        Args:
            values (list of str): The enum values that should be selected.
        """
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
        """Return what values of the enumeration have been selected.
        
        Return (list of str): A comma seperated list of the selected values.
        """
        return ",".join(sorted(self._values))


class PypeItCustomEditorDelegate(QStyledItemDelegate):
    """Custom item delegate for rows in a PypeItMetadataView."""
    def __init__(self, parent):
        self.metadata_view=parent
        super().__init__(parent)


    def editorClosed(self,editor, accepted):
        """Signal handler that is notified when the PypeItEnumListEditor is closed.
        
        Args:
            editor (PypeItEnumListEditor): The editor that was closed.
            accepted (bool): True if the value was accepted, false if it was canceled.
        """

        # Notify any clients to accept or revert any cached values. This signal is inherited from
        # parent classes.
        if accepted:
            self.closeEditor.emit(editor, QAbstractItemDelegate.EndEditHint.SubmitModelCache)
            self.setModelData(editor, editor.index.model(), editor.index)
        else:
            self.closeEditor.emit(editor, QAbstractItemDelegate.EndEditHint.RevertModelCache)

        
    def paint(self, painter, option, index):
        """
        Overridden version of paint for painting items in the PypeItMetadataView.
        """
        # This method is overridden because we want to display commented out items
        # as if the are disabled, without actually disabling them. This allows them to
        # be selected and uncommented out.
        
        # Get the source model of ProxyModel PypeItMetadataView uses to sort items.
        # This source model is a PypeItMetadataModel
        model=index.model().sourceModel()
        source_index = index.model().mapToSource(index)
        if model.isCommentedOut(source_index):
            # Paint as disabled
            option.state &= ~QStyle.StateFlag.State_Enabled

        super().paint(painter,option, index)

    def createEditor(self, parent,  option, index):
        """
        Creates an editor widget for an item in the metadata table. This will be a
        PypeItEnumListEditor for the columns that use one, or the Qt default for other editable columns.
        Overriden from QStyledItemDelegate.

        Args:
            parent (QWidget): The parent widget of the new editor.
            option (QtWidgets.QStyleOptionViewItem): Additional options for the editor.
            index (QModelIndex): The index of the table cell being edited.
        """
        model = index.model().sourceModel()
        
        column_name = model.getColumnNameFromNum(index)

        if column_name == "frametype":
            msgs.info("Creating enum list editor for frametype")
            editor= PypeItEnumListEditor(parent=parent, index=index, num_lines=5, allowed_values=model.getAllFrameTypes())
            editor.closed.connect(self.editorClosed)
            return editor
        
        msgs.info(f"Creating default editor for {column_name}")
        return super().createEditor(parent, option, index)
    
    def setEditorData(self, editor, index):
        """Sets the data being edited in the editor. Overriden from QStyledItemDelegate.

        Args:
            editor (QWidget):    The editor widget (created by createEditor)
            index (QModelIndex): The index of the item being edited.
        """
        if isinstance(editor, PypeItEnumListEditor):
            msgs.info(f"Setting editor data {index.data(Qt.EditRole)}")
            editor.setSelectedValues(index.data(Qt.EditRole))
        else:
            msgs.info("Setting default editor data")
            super().setEditorData(editor, index)

    def setModelData(self,editor,model,index):
        """Sets the edited data in the model post editing. Overriden from QStyledItemDelegate.
        
        Args:
            editor (QWidget):           The editor widget (created by createEditor).
            model (QAbstractItemModel): The model being edited.
            index (QModelIndex):        The index of the item being edited.
        """
        if isinstance(editor, PypeItEnumListEditor):
            msgs.info(f"Setting choice model data: {editor.selectedValues()}")
            model.setData(index, editor.selectedValues())
        else:
            msgs.info("Setting default model data")
            super().setModelData(editor,model,index)

    def updateEditorGeometry(self, editor, option, index):
        """Sets the editor's position and size in the GUI. Overriden from QStyledItemDelegate.
        
        Args:
            editor (QWidget):                        The editor widget (created by crateEditor). This widgets geometry
                                                     is set by this method.
            model (QAbstractItemModel):              The model being edited
            option (QtWidgets.QStyleOptionViewItem): Options object containing the recommended rectangle for the editor.
            index (QModelIndex):                     The index of the item being edited.
        """
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

            # Because the parent may be in a scrollable, the x,y could be negative,
            # or underneath the header row. Set those so that the editor doesn't appear outside of the widget
            if editor_x < 0:
                editor_x = 0

            min_y = self.metadata_view.verticalHeader().sectionSize(0)
            if editor_y < min_y:
                editor_y = 0

            # Adjust the editor'x upper left corner so that the editor is vislbe for
            # cells along the bottom or right of the parent table
            right_x  = self.metadata_view.viewport().geometry().bottomRight().x() - self.metadata_view.viewportMargins().right()
            bottom_y = self.metadata_view.viewport().geometry().bottomRight().y() - self.metadata_view.viewportMargins().bottom()

            # The bottom x,y of the view port is measured without the margins, but the editor is placed relative to those
            # margins, so we have to include the left/top margins in the below calculations
            if editor_x + self.metadata_view.viewportMargins().left() + editor_width > right_x:
                editor_x = right_x - (editor_width + self.metadata_view.viewportMargins().left())

            if editor_y + self.metadata_view.viewportMargins().top() + editor_min_size.height() > bottom_y:
                editor_y = bottom_y - (editor_min_size.height() + self.metadata_view.viewportMargins().top())

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

        # Set to a minimum number of rows high so the frame type editor has enough space
        if model.rowCount() > 0:
            row_height = self.verticalHeader().sectionSize(0)
        else:
            row_height = self.verticalHeader().defaultSectionSize()

        # We use 11 rows, 1 for the header, and 10 data rows. This seems to give an adaquate buffer to the frame type editor.
        min_height = (self.contentsMargins().top() + self.contentsMargins().bottom() + 
                      self.horizontalScrollBar().sizeHint().height() +
                      11*row_height)
        msgs.info(f"current min_height/height/hint h: {self.minimumHeight()}/{self.height()}/{self.sizeHint().height()}, scrollbar hint h {self.horizontalScrollBar().sizeHint().height()}, currentmargin top/bottom: {self.contentsMargins().top()}/{self.contentsMargins().bottom()} hdr min_height/height/hint h: {self.horizontalHeader().minimumHeight()}/{self.horizontalHeader().height()}/{self.horizontalHeader().sizeHint().height()}")
        msgs.info(f"rowHeight: {row_height} current min_height {self.minimumHeight()} new min_height {min_height}")
        if min_height > self.minimumHeight():
            self.setMinimumHeight(min_height)

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
        config (dict):     The name/value pairs for the configuration keys defined by the spectrograph.
        lines_to_display (int):     How many lines to display before scrolling.
        parent (QWidget, Optional): The parent widget, defaults to None.
    """
    def __init__(self, spec_name, config, lines_to_display, parent=None):

        super().__init__(parent=parent)        
        self.setTitle(self.tr("Setup"))

        # Set our configuration, adding the spectrograph name as the first item.
        self.lines_to_display = lines_to_display

        # Put everything in a scroll area
        layout = QHBoxLayout(self)
        self._scroll_area = QScrollArea(parent=self)
        self._scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        # A widget using a form layout to display the spectrograph + configuration values
        self._form_widget = QWidget()
        self._form_widget_layout = QFormLayout(self._form_widget)

        # Add margins to the right to avoid needing a horizontal scroll bar
        fm = self.fontMetrics()
        max_char_width = fm.maxWidth()
        self._form_widget_layout.setContentsMargins(0, 0, max_char_width, 0)

        # Keep this from expanding too large.
        self._form_widget_layout.setSizeConstraint(QLayout.SizeConstraint.SetMinimumSize)
        
        # Add a "Spectrograph" line to the form layout
        self._config_labels = {}
        label = QLabel(spec_name)
        self._config_labels['Spectrograph'] = label
        self._form_widget_layout.addRow('Spectrograph', label)

        # Add values to form layout
        for key, value in config.items():
            label = QLabel(str(value))
            self._form_widget_layout.addRow(key, label)
            self._config_labels[key] = label

        self._scroll_area.setWidget(self._form_widget)

        # Set the minimum width of the form widget
        self._form_widget.setMinimumWidth(self._getMinWidth())


        # Figure out the correct height for this panel, so that only the spectrograph and self.number_of_lines
        # config keys are visible

        # Find the minimum height of the form widget needed to hold the number of lines to display
        msgs.info(f"font height: {fm.height()} vertical spacing {self._form_widget_layout.verticalSpacing()}")
        min_fw_height = self._form_widget_layout.verticalSpacing()*(lines_to_display-1) + fm.height()*lines_to_display

        # The height of this panel is that height plus the margins + the group box title
        scroll_area_margins = self._scroll_area.contentsMargins()
        group_box_margins = self.contentsMargins()
        form_widget_margins = self._form_widget.contentsMargins()
        self.setFixedHeight(min_fw_height + 
                            fm.height()   +  # Group Box Title
                            group_box_margins.top()   + group_box_margins.bottom() +
                            scroll_area_margins.top() + scroll_area_margins.bottom() +
                            form_widget_margins.top() + form_widget_margins.bottom())

        # Set to fixed sizing policy
        policy = QSizePolicy()
        policy.setHorizontalPolicy(QSizePolicy.Minimum)
        policy.setVerticalPolicy(QSizePolicy.Fixed)
        policy.setControlType(QSizePolicy.DefaultType)
        self.setSizePolicy(policy)
        
        layout.addWidget(self._scroll_area)

    def setNewValues(self, config_dict: dict) -> None:
        """Update the panel to display new configuration values.
        
        Args:
            config_dict: A dict of the new values.
        """
        for key, value in config_dict.items():
            label = self._config_labels.get(key,None)
            if label is None:
                label = QLabel(str(value))
                self._config_labels[key] = label
                self._form_widget_layout.addRow(key,label)
            else:
                label.setText(str(value))

        # Reset the minimum width for the new values
        self._form_widget.setMinimumWidth(self._getMinWidth())


    def _getMinWidth(self) -> int:
        """Calculate the minimum width needed to display the configuration values."""

        # We calculate our width to just fit around the text so we can stop the scroll area from
        # covering it with a scrollbar.
        fm = self.fontMetrics()
        max_key_width = 0
        max_value_width = 0
        for key in self._config_labels:
            value = self._config_labels[key].text()
            key_width = fm.horizontalAdvance(key)
            value_width = fm.horizontalAdvance(value)

            if key_width > max_key_width:
                max_key_width = key_width

            if value_width > max_value_width:
                max_value_width = value_width

        margins = self._form_widget_layout.contentsMargins()

        min_width = self._form_widget_layout.horizontalSpacing() + max_key_width + max_value_width + margins.left() + margins.right()

        # Account for the scroll bar if needed
        if len(self._config_labels) > self.lines_to_display:
            if self._scroll_area.verticalScrollBar():
                min_width += self._scroll_area.verticalScrollBar().sizeHint().width()
        msgs.info(f"new minWidth: {min_width} max key width: {max_key_width} max_value width {max_value_width} horizontal spacing {self._form_widget_layout.horizontalSpacing()} margins left: {margins.left()} margins right: {margins.right()}")
        return min_width

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
        """str: The name of the tab."""
        return self._name

    @property
    def state(self):
        """ModelState: The state of the tab."""
        return self._state

    @property
    def closeable(self):
        """bool: Whether the tab can be closed."""
        return self._closeable

class PypeItFileView(TabManagerBaseTab):
    """Widget for displaying the information needed for one pypeit file. This includes
    the spectrograph, the configuration keys and values, the files that share that 
    configuration, and the PypeIt parameters.

    Args:
        (pypeit.setup_gui.model.PypeItFileModel): The model representing all the information needed for a .pypeit file.
        (pypeit.setup_gui.model.PypeItFileController): The controller for managing the user's interaction with a PypeIt file)
    """


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
        self.config_panel = ConfigValuesPanel(model.spec_name, model.config_values, 5, parent=self)
        third_row_layout.addWidget(self.config_panel)

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

        # Create a group box and table view for the file metadata table
        file_group = QGroupBox(self.tr("File Metadata"))
        file_group_layout = QVBoxLayout()
        self.file_metadata_table = PypeItMetadataView(self, model.metadata_model, controller.getMetadataController(model.metadata_model))
        file_group_layout.addWidget(self.file_metadata_table)
        file_group.setLayout(file_group_layout)        
        layout.addWidget(file_group)

        
        # Stretch the metadata and params rows more than the filename and config_key rows
        layout.setStretch(2,4)
        layout.setStretch(3,10)
        layout.setStretch(4,10)

        self.model.stateChanged.connect(self.update_from_model)


    def update_from_model(self):
        """
        Signal handler that updates view when the underlying model changes.
        """
        # update the filename if it changed from saving
        self.filename_value.setText(self.model.filename)

        # Update the config values
        self.config_panel.setNewValues(self.model.config_values)

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
        controller (:class:`pypeit.setup_gui.controller.PypeItObsLogController`): Controller object for the PypeIt Setup GUI.
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
        self.spectrograph.addItems(available_spectrographs)
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
        """Signal handler that removes raw data paths from the obslog"""
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
            if str_input.lower() in  available_spectrographs:
                return QValidator.Acceptable, str_input, int_input
            else:
                for spectrograph in  available_spectrographs:
                    if spectrograph.startswith(str_input.lower()):
                        return QValidator.Intermediate, str_input, int_input
        return QValidator.Invalid, str_input, int_input

class TabManagerWidget(QTabWidget):
    """
    Widget which manages tabs for the obslog and pypeit files.
    It extends the QTabWidget functionality by allowing tabs
    to be added or removed and by displaying the tab's name with a 
    "*" when a tab hasn't been saved.

    Args:
        model (:class:`pypeit.setup_gui.model.SetupGUIStateModel`): Model object for a PypeIt Setup GUI.
        controller (:class:`pypeit.setup_gui.model.SetupGUIControllerl`): Controller object for the PypeIt Setup GUI.
    """

    tabCreateRequest = Signal()
    """Signal sent when the wants to add a new tab."""

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
        """Signal handler that detects a click on the "+" tab widget and sends that as a tabCreateRequest signal."""
        if index == self.count() - 1:
            # Create as new tab. The new tab model will send the signals needed
            # to create the view and add it
            self.tabCreateRequest.emit()

    def checkIfNewTabCurrent(self, index):
        """Signal handler thast prevents the "+" tab, from being the current tab."""
        # Try to prevent the + tab from being visible
        if self.count() > 1:
            if self.widget(index).name == "+":
                # Make the first tab visible instead
                self.setCurrentIndex(0)

    def addNewTab(self, tab):
        """ Insert a new before the "+" tab. If this is the first tab being entered
        (ie the + tab itself), it is appended
        
        Args:
            tab (TabManagerBaseTab): The new tab widget to add.

        Rerturns:
            int: The index of the newly inserted tab.
        """
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
        """Close the tab with the given name.
        
        Args:
            tab_name (str): The name of the tab to close.
        """
        try:
            index = self._tabNames.index(tab_name)
        except ValueError :
            msgs.warn(f"Failed to find tab named {tab_name} in list.")
            return
        tab = self.widget(index)
        if tab.closeable:
            tab.stateChanged.disconnect(self.updateTabText)
        self.removeTab(index)
        del self._tabNames[index]

    def updateTabText(self, tab_name, tab_state):
        """Update a tab's text when it state changes.
        
        Args:
            tab_name (str): The name of the tab to update.
            tab_state (ModelState): The model state of the data in the tab. If this is NEW or UNCHANGED,
                                    the tab's text will be  "*" + tab_name. Otherwise it will be
                                    tab_name.
        """
        try:
            index = self._tabNames.index(tab_name)
        except ValueError :
            msgs.warn(f"Failed to find tab named {tab_name} in list.")
            return

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
        self.tab_widget.tabCloseRequested.connect(self._closeRequest)

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

    def _closeRequest(self, index):
        """Called when the user tries to close a tab"""
        tab = self.tab_widget.widget(index)
        if tab.closeable:
            self.controller.close(tab.model)

    def _showLog(self):
        """Signal handler that opens the log window."""
        if self._logWindow is not None:
            self._logWindow.activateWindow()
            self._logWindow.raise_()
        else:
            self._logWindow = LogWindow(self.model.log_buffer)
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

