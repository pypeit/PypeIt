"""
Utilities for creating and displaying dialogs in Qt.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

# Allow runtime evaluation of type annotations, specifically to allow class methods to return their own class
from __future__ import annotations

import enum
from typing import Optional,Union
from pathlib import Path
from dataclasses import dataclass

from pypeit import msgs

from qtpy.QtWidgets import QFileDialog, QMessageBox, QCheckBox, QDialog, QWidget
from qtpy.QtCore import QStringListModel, QSettings

@dataclass
class FileType:
    name:str
    extension:str


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


class FileDialog:
    """Opens a file dialog for either opening or saving files.  This encapsulates
    the logic to use QFileDialog in a way that can be used by the view or controller and
    can be easilly patched by unit tests.
    
    Args:
        parent:      Parent of the dialog widget.
        caption:     Caption for the dialog.
        file_mode:   The file mode for the dialog.
        filter:      The filter to use for files in the file dialog. For example: "PypeIt input files (\*.pypeit)"
        history:     The list of paths to populate the history of the dialog. New paths selected are added to this history.
        save:        Whether the dialog is a save dialog. If False it is treated as an open dialog. Defaults to False.
        ask_for_all: Whether the dialog should present a "Use this location for everything" option
                                           when saving. Defaults to False.
     """
    def __init__(self, 
                 parent : QWidget, 
                 caption : str, 
                 file_mode : QFileDialog.FileMode, 
                 file_type : Optional[FileType] =None,
                 default_file : Optional[Union[str,Path]] = None,
                 history : Optional[QStringListModel] = None, 
                 save : bool = False, 
                 ask_for_all : bool = False):
        
        if file_type is not None:
            filter = f"{file_type.name} (*{file_type.extension})"
        else:
            filter = None
        self._dialog = QFileDialog(parent, caption=parent.tr(caption), filter=filter)
        self._ask_for_all = ask_for_all

        self._dialog.setOption(QFileDialog.Option.DontUseNativeDialog)
        if save:
            self._dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
        else:
            self._dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptOpen)

        self._dialog.setFileMode(file_mode)


        self._history = history
        if history is not None and history.rowCount() > 0:
            self._dialog.setDirectory(history.stringList()[-1])
            self._dialog.setHistory(history.stringList())


        if default_file is not None:
            # Set the starting directory to that of the default file, if it has a directory
            if isinstance(default_file,str):
                default_file = Path(default_file)
            starting_dir = default_file.parent
            if starting_dir.parent.name != "":
                self._dialog.setDirectory(str(starting_dir))
            self._dialog.selectFile(str(default_file))

        # Adding a checkbox like this only works for non-native dialogs. And if the 
        # default layout of the QFileDialog ever changes this could stop working.
        if self._ask_for_all:
            sd_layout = self._dialog.layout()
            row_count = sd_layout.rowCount()
            self.use_for_all_checkbox = QCheckBox(self._dialog.tr("Use this location for everthing."), self._dialog)
            sd_layout.addWidget(self.use_for_all_checkbox, row_count, 0)

        self.selected_path = None

    def show(self) -> DialogResponses:
        """Show a modal file dialog. If the dialog is accepted, the results will be saved
        in the selected_path and response attributes.
        
        Returns: 
            DialogResponses: CANCEL, ACCEPT_FOR_ALL, or ACCEPT.
        """
        result = self._dialog.exec()
        if result ==  QDialog.Rejected:
            self.response = DialogResponses.CANCEL
            self.selected_path = None
            return DialogResponses.CANCEL

        path = Path(self._dialog.selectedFiles()[0])
        self.selected_path = str(path)
        if self._history is not None:
            # Append new selections to history
            # But we only want the path, not the full filename
            if not path.is_dir():
                path = path.parent
            if str(path) not in self._history.stringList():
                row_index = self._history.rowCount()
                self._history.insertRows(row_index, 1)
                self._history.setData(self._history.index(row_index, 0),str(path))

        if self._ask_for_all and self.use_for_all_checkbox.isChecked():
            return DialogResponses.ACCEPT_FOR_ALL

        return DialogResponses.ACCEPT
            
    @classmethod
    def create_open_file_dialog(cls, parent : QWidget, caption : str, file_type : FileType, history_group : str = "OpenFile") -> FileDialog:
        """Creates a dialog to prompt the user for an existing file.
        
        Args:
            parent: The parent widget of the pop-up dialog
            caption: A caption for the dialog
            filter:  A QFileDialog filter for a file to open. For example: "Text Files (`*`.txt)"
            history_group: The group in the applications persistent settings to persist the history.
                           Defaults to "OpenFile"
                               

        Returns:
            FileDialog: The created file dialog.
        """
        # Get history from settings
        history = PersistentStringListModel(history_group, "History")

        return cls(parent, caption, QFileDialog.ExistingFile, file_type=file_type, history=history)

    @classmethod
    def create_save_location_dialog(cls, parent : QWidget, config_name : str, prompt_for_all : bool =False, history_group="OutputDirectory") -> FileDialog:
        """Creates a dialog requesting the user select a location
        to save a file. A history is maintained.

        Args:
            parent:         The parent widget of the pop-up dialog
            config_name:    The name of the configuration being saved to a pypeit file.
            prompt_for_all: Whether to prompt the user if this save location should
                            apply to all subsequent files in the operation.
            history_group:  The group in the applications persistent settings to persist the history.
                            Defaults to "OutputDirectory"
                               

        Returns:
            The FileDialog object used to prompt the user. It's response attribute will be either CANCEL, ACCEPT, or ACCEPT_FOR_ALL.
            It's selected_path members will be set to the chosen path or None of the dialog was canceled.
        """

        # Get history
        history = PersistentStringListModel(history_group, "History")

        # Create the dialog.
        return cls(parent, parent.tr(f"Select location to save tab {config_name}."),
                   QFileDialog.Directory, history=history, save=True, ask_for_all=prompt_for_all)


def display_error(parent : QWidget, message: str) -> None:
    """
    Display an error message pop up dialog to the user.

    Args:
        parent: The parent widget of the pop-up dialog
        message: The message to display.
    """
    msgs.warn(message) # Make sure the message also goes to the logs
    QMessageBox.warning(parent, "PypeIt Setup Error", message, QMessageBox.Ok)

def prompt_to_save(parent : QWidget) -> DialogResponses:
    """Prompt the user if they want to save any unsaved data.

    Args:
        parent: The parent widget of the pop-up dialog

    Returns:
        DialogResponses: 
            ``SAVE`` if the user wants to save, ``ACCEPT`` if they want to continue 
            without saving, ``CANCEL`` if they want to cancel the current operation 
            and to not lose any data.
    """
    response = QMessageBox.warning(parent, parent.tr("PypeIt Setup"), parent.tr("There are unsaved changes to PypeIt setup files.\nDo you want to save then?"),
                                   QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel, QMessageBox.Save)
    if response == QMessageBox.Save:
        return DialogResponses.SAVE
    elif response == QMessageBox.Discard:
        return DialogResponses.ACCEPT
    else:
        return DialogResponses.CANCEL

