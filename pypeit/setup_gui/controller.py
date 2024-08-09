"""
The controller portion of the PypeIt Setup GUI. The classes in this module are responsible for
acting on user input, running background tasks, and returning information to the user.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import traceback
import sys
import threading
from datetime import datetime
import re
import io
from pathlib import Path
from functools import partial
from contextlib import contextmanager
from qtpy.QtCore import QCoreApplication, Signal, QMutex, QTimer

# TODO: datetime.UTC is not defined in python 3.10.  Remove this when we decide
# to no longer support it.
try:
    __UTC__ = datetime.UTC
except AttributeError as e:
    from datetime import timezone
    __UTC__ = timezone.utc

from qtpy.QtCore import QObject, Qt, QThread
from qtpy.QtGui import QKeySequence
from qtpy.QtWidgets import QAction

from pypeit.setup_gui.view import SetupGUIMainWindow, PypeItFileView, DialogResponses
from pypeit.setup_gui.text_viewer import TextViewerWindow
from pypeit.setup_gui.dialog_helpers import prompt_to_save, display_error, FileDialog, FileType
from pypeit.setup_gui.model import PypeItSetupGUIModel, ModelState, PypeItFileModel
from pypeit import msgs

from pypeit.display import display
from pypeit import io as pypeit_io


@contextmanager
def lock_qt_mutex(mutex):
    """Context manager to allow locking Qt :class:`QMutex` objects using 'with'."""
    mutex.lock()
    try:
        yield mutex
    finally:
        mutex.unlock()

class OpCanceledError(Exception):
    """Exception thrown when a background operation has been canceled."""
    def __init__(self):
        super().__init__()

class OperationThread(QThread):
    """Thread to run a background operation."""

    completed = Signal(bool, tuple)
    """Signal send the operation has completed."""


    def __init__(self):
        super().__init__()
        self._operation = None
        self._max_progress = None
        self._mutex = QMutex()

    def run(self):
        """Runs an operation in a background thread."""
        canceled = False
        exc_info = (None, None, None)
        try:
            self._operation.run()            
        except OpCanceledError:
            canceled=True
        except Exception:
            exc_info=sys.exc_info()

        self.completed.emit(canceled, exc_info)        

    def _cancel_op(self):
        """Cancels an in progress background operation when the user cancels the progress dialog."""
        self.requestInterruption()

    def _op_progress(self, max_progress, progress_message=None):
        """Signal handler that receives progress information from a background operation
        as it proceeds. It passes this to the view to increase the value showing in the
        progress dialog."""

        with lock_qt_mutex(self._mutex):
            create_progress = False
            mp = None
            if self._operation is not None:
                if self._max_progress is None:
                    if max_progress is not None:
                        self._max_progress = max_progress
                        create_progress = True
                mp = self._max_progress
        
        if create_progress:
            # If we've just initialized the max progress, create the progress dialog
            SetupGUIController.main_window.create_progress_dialog(self._operation.name, max_progress, self._cancel_op)
            
        # Ignore the progress if there's no max progress yet
        if mp is not None:            
            SetupGUIController.main_window.show_operation_progress(increase=1, message = progress_message)

    def _op_complete(self, canceled, exc_info):
        """Signal handler that is notified when a background operation completes.
        
        Args:
            canceled (bool): Whether or not the operation was canceled.
            exc_info (tuple): The exception information if the operation failed. None if it succeeded
        """
        msgs.info("Op complete")
        with lock_qt_mutex(self._mutex):
            if self._operation is not None:
                operation = self._operation
                self._operation = None
                self._max_progress = None
            else:
                operation = None

        if operation is not None:
            operation.progressMade.disconnect(self._op_progress)
            SetupGUIController.main_window.operation_complete()
            operation.postRun(canceled, exc_info)            
    

    def startOperation(self, operation):
        """
        Start a background operation.
        
        Args:
            operation (MetadataOperation): The MetadataOperation to start in the background thread.
        """
        self._operation = operation
        if operation.preRun():
            operation.progressMade.connect(self._op_progress, type=Qt.QueuedConnection)
            self.completed.connect(self._op_complete, type=Qt.QueuedConnection)
            self.start()

class MetadataOperation(QObject):

    """Base class for Metadata operations that take long enough that they should take place in a background thread.
    
    Args:
        model (PypeItSetupGUIModel): The PypeIt Setup GUI's model object.
    """

    progressMade = Signal(int, str)
    """Signal emitted emit when progress has been made. This will be reflected in the view's progress dialog."""


    def __init__(self, name, model):
        super().__init__()
        self._model=model
        self.name = name
        self._max_progress = None

    def preRun(self):
        """
        Perform setup required before running the operation. This involves watching the log
        for files being added to the metadata.
        """
        building_metadata_re = re.compile(r"Building metadata for (\d+) ")
        self._model.log_buffer.watch("building_metadata", building_metadata_re, self._buildingMetadata)
        
        added_metadata_re = re.compile(r"Adding metadata for (.*)$")
        self._model.log_buffer.watch("added_metadata", added_metadata_re, self._addedMetadata)
        self._model.closeAllFiles()
        return True

    def _buildingMetadata(self, name, match):
        """Callback used to find the total number of files being read when building metadata."""
        self._max_progress = int(match.group(1))
        msgs.info(f"Found max progress {self._max_progress}")

    def _addedMetadata(self, name, match):
        """Callback used to report progress on reading files when building metadata."""
        if QThread.currentThread().isInterruptionRequested():
            raise OpCanceledError()

        self.progressMade.emit(self._max_progress, match.group(1))


    def postRun(self, canceled, exc_info):
        """Clean up steps after the operations has run.
        
        Args:
            canceled (bool):  True if the operation was canceled. 
            exc_info (tuple): The exception information (as returned by sys.exc_info()) for any errors that occurred.
        """
        self._model.log_buffer.unwatch("added_metadata")
        self._model.log_buffer.unwatch("building_metadata")

        if exc_info[0] is not None:
            traceback_string = "".join(traceback.format_exception(*exc_info))
            msgs.warn(f"Failed to {self.name.lower()}:\n" + traceback_string)
            display_error(SetupGUIController.main_window, f"Failed to {self.name.lower()} {exc_info[0]}: {exc_info[1]}")
            self._model.reset()
        elif canceled:
            self._model.reset()

    def run(self):
        """Performs the steps of the MetadataOperation. This is an abstract method overridden by child classes."""
        pass

class SetupOperation(MetadataOperation):
    """ Background operation to run pypeit setup.

        Args:
            model (PypeItSetupGUIModel): The PypeIt Setup GUI's model object.
    """
    def __init__(self, model):
        super().__init__("Run Setup", model)

    def run(self):
        """
        Runs pypeit_setup on the current raw data directories.
        """
        self._model.run_setup()

class OpenFileOperation(MetadataOperation):
    """
    Background operation to open a PypeIt file

    Args:
        model (PypeItSetupGUIModel):
            The PypeIt Setup GUI's model object.
        file (): 
            The file to open.
    """

    def __init__(self, model, file):
        super().__init__("Open File", model)
        self._file = file

    def run(self):
        """
        Opens a pypeit file and reads metadata for all of the files in it.
        """
        self._model.open_pypeit_file(self._file)

class MetadataReadOnlyAction(QAction):
    """An action on a PypeItMetadataModel that is read only and therefore can be performed on the ObsLog.
    These actions can be triggered by a button, menu, option or keyboard shortcut.
    
    Arguments:
        controller(PypeItMetadataController): 
            The controller for the PypeItMetadataModel/PypeItMetadataView MVC triplet.
        menu_text(str):
            The text name for the menu/button that triggers the action.
        handler (:obj:`collections.abc.Callable`):
            The signal handler to enact the action. This receives the "triggered" event
            from the parent class.
        shortcut (QtGui.QKeySequence.StandardKey, Optional):
            The keyboard shortcut to initiate the action.
    """
    def __init__(self, controller, menu_text, handler, shortcut=None):
        super().__init__(menu_text)
        self.triggered.connect(handler)
        if shortcut is not None:
            self.setShortcut(shortcut)
        self._controller=controller

    def updateEnabledStatus(self):
        """Enable/disable the action based on whether any metadata rows are selected."""
        self.setEnabled(self._controller._view is not None and len(self._controller._view.selectedRows()) > 0)

class MetadataWriteAction(QAction):
    """An action on a PypeItMetadataModel that can change the file metadata and therefore can only be 
    performed on a PypeItFileModel.
    These actions can be triggered by a button, menu, option or keyboard shortcut.

    Arguments:
        controller(PypeItMetadataController): 
            The controller for the PypeItMetadataModel/PypeItMetadataView MVC triplet.
        menu_text(str):
            The text name for the menu/button that triggers the action.
        handler (:obj:`collections.abc.Callable`):
            The signal handler to enact the action. This receives the "triggered" event
            from the parent class.
        shortcut (QtGui.QKeySequence.StandardKey, Optional):
            The keyboard shortcut to initiate the action.
    """

    def __init__(self, controller, menu_text, handler, shortcut=None):
        super().__init__(menu_text)
        self.triggered.connect(handler)
        if shortcut is not None:
            self.setShortcut(shortcut)
        self._controller=controller

    def updateEnabledStatus(self):
        """Enable/disable the action based on whether any metadata rows are selected."""

        self.setEnabled(self._controller._is_pypeit_file and
                        self._controller._view is not None and 
                        len(self._controller._view.selectedRows()) > 0)

class MetadataPasteAction(QAction):
    """An action for pasting metadata into a PypeItMetadataModel object.
    
    Arguments:
        controller(PypeItMetadataController): 
            The controller for the PypeItMetadataModel/PypeItMetadataView MVC triplet.
        menu_text(str):
            The text name for the menu/button that triggers the action.
        handler (:obj:`collections.abc.Callable`):
            The signal handler to enact the action. This receives the "triggered" event
            from the parent class.
        shortcut (QtGui.QKeySequence.StandardKey, Optional):
            The keyboard shortcut to initiate the action.
    """

    def __init__(self, controller, menu_text, handler, shortcut=None):
        super().__init__(menu_text)
        self.triggered.connect(handler)
        if shortcut is not None:
            self.setShortcut(shortcut)
        self._controller=controller

    def updateEnabledStatus(self):
        """Enable/disable the action based on whether any metadata rows are selected AND
        there are rows to paste in the clipboard."""

        if self._controller._is_pypeit_file:
            if SetupGUIController.model.clipboard.rowCount() > 0:
                self.setEnabled(True)
            else:
                self.setEnabled(False)
        else:
            self.setEnabled(False)

class PypeItMetadataController(QObject):
    """A Controller object for performing actions iniitiated by the user file metadata.
    Part of a MVC triplet involving PypeItMetadataModel/PypeItMetadataController/PypeItMetadataView.

    Args:
        model (:obj:`pypeit.setup_gui.model.PypeItMetatadataModel`):
            The model this controller acts with.

        is_pypeit_file (bool): 
            True if the model is for a PypeItFileModel (that is writeable model), False if it is 
            from a PypeItObsLog model (read only)         
    """
    def __init__(self, model, is_pypeit_file):
        super().__init__()
        self._model = model
        self._view = None
        self._is_pypeit_file = is_pypeit_file
        self._windows = {}
        self.next_window_id = 1

        # Build actions
        self._action_list = [MetadataReadOnlyAction(self, "View File",   self.view_file),
                             MetadataReadOnlyAction(self, "View Header", self.view_header),
                             MetadataReadOnlyAction(self, "Copy",        self.copy_metadata_rows,       shortcut=QKeySequence.StandardKey.Copy),
                                MetadataWriteAction(self, "Cut",         self.cut_metadata_rows,        shortcut=QKeySequence.StandardKey.Cut),
                                MetadataPasteAction(self, "Paste",       self.paste_metadata_rows,      shortcut=QKeySequence.StandardKey.Paste),
                                MetadataWriteAction(self, "Comment Out", self.comment_out_metadata_rows),
                                MetadataWriteAction(self, "Uncomment",   self.uncomment_metadata_rows),
                                MetadataWriteAction(self, "Delete",      self.delete_metadata_rows,     shortcut=QKeySequence.StandardKey.Delete) ]
        SetupGUIController.model.clipboard.modelReset.connect(self.updatedEnabledActions)
        self.updatedEnabledActions()

    def getActions(self, parent):
        """Returns the actions that this controller supports.
        
        Returns: 
            list of QAction: List of the actions that can be performed on the PypeItMetadataModel.
        """
        return self._action_list
            
    def setView(self, view):
        """Set the view that is responsible for displaying and receiving input for
        the PypeItMetadataModel.
        
        Args:
            view (PypeItMetadataView): The view.
        """
        self._view = view
        view.selectionUpdated.connect(self.updatedEnabledActions)

    def updatedEnabledActions(self):
        """Updates which actions are enabled/disabled."""
        for action in self._action_list:
            action.updateEnabledStatus()
                    
    def view_file(self):
        """View the selected files in the metadata using Ginga."""
        row_indices = self._view.selectedRows()
        if len(row_indices) > 0:

            # Make sure ginga is available
            try:
                display.connect_to_ginga(raise_err=True, allow_new=True)
            except Exception as e:
                display_error(SetupGUIController.main_window, f"Could not start ginga to view FITS files: {e}")
                msgs.warn(f"Failed to connect to ginga:\n" + traceback.format_exc())

            
            # Display each file in its own ginga tab
            for indx in row_indices:
                metadata_row = self._model.metadata[indx]
                # Make sure to strip comments off commented out files
                file = Path(metadata_row['directory'], metadata_row['filename'].lstrip('# '))
                if not file.exists():
                    display_error(SetupGUIController.main_window, f"Could not find {file.name} in {file.parent}.")
                    return

                try:
                    img = self._model.spectrograph.get_rawimage(str(file), 1)[1]
                except Exception as e:
                    display_error(SetupGUIController.main_window, f"Failed to read image {file.name}: {e}")
                    msgs.warn(f"Failed get raw image:\n" + traceback.format_exc())

                try:
                    display.show_image(img, chname = f"{file.name}")
                except Exception as e:
                    display_error(SetupGUIController.main_window, f"Failed to send image {file.name} to ginga: {e}")
                    msgs.warn(f"Failed send image to ginga:\n" + traceback.format_exc())

    def view_header(self):
        """ Display the header of one or more selected files in the metadata.
        """
        row_indices = self._view.selectedRows()
        if len(row_indices) > 0:            
            # Display each file in its window
            for indx in row_indices:
                metadata_row = self._model.metadata[indx]
                # Make sure to strip comments off commented out files
                file = Path(metadata_row['directory'], metadata_row['filename'].strip('# '))
                header_string_buffer = io.StringIO()
                try:
                    with pypeit_io.fits_open(file) as hdul:
                        for i, hdu in enumerate(hdul):
                            print(f"\n\n# HDU {i} Header from {file.name}\n",file=header_string_buffer)
                            hdu.header.totextfile(header_string_buffer)
                except Exception as e:
                    display_error(SetupGUIController.main_window, f"Failed to read header from file {file.name} in {file.parent}: {e}")
                    msgs.warn(f"Failed to read header from {file}:\n" + traceback.format_exc())
                    return
                header_string_buffer.seek(0)
                window = TextViewerWindow(title=f"{file.name} Header", width=80, height=50,start_at_top=True, filename=file.parent / (file.name+".txt"), text_stream=header_string_buffer)
                self._windows[self.next_window_id] = window
                window.closed.connect(partial(self._window_closed, id=self.next_window_id))
                self.next_window_id+=1
                window.show()

    def _window_closed(self, id : int) -> None:
        """Clean up when a header viewer window is closed"""
        del self._windows[id]

    def copy_metadata_rows(self) -> bool:
        """Copy metadata rows into the clipboard.
        
        Return: 
            True if rows were copied, False if there were no rows to copy
        """
        if self._view is None:
            return False

        row_indices = self._view.selectedRows()
        msgs.info(f"Copying {len(row_indices)} rows to the clipboard.")
        if len(row_indices) > 0:
            row_model = self._model.createCopyForRows(row_indices)
            SetupGUIController.model.clipboard = row_model
            return True
        return False

    def cut_metadata_rows(self) -> bool:
        """Move metadata rows from the PypeItMetadataModel to the clipboard.
        
        Return: 
            True if rows were removed from the metadata and copied into the
            clipboard. False if no rows were copied or removed.
        """
        if self.copy_metadata_rows():
            return self.delete_metadata_rows()
        return False


    def paste_metadata_rows(self):
        """Insert rows from the clipboard into the PypeItMetadataModel"""
        clipboard = SetupGUIController.model.clipboard
        if clipboard.rowCount() > 0:
            try:
                msgs.info(f"Pasting {clipboard.rowCount()} rows")
                self._model.pasteFrom(clipboard)
            except Exception as e:
                traceback_string = "".join(traceback.format_exc())
                msgs.warn(f"Failed to paste metadata rows:\n" + traceback_string)
                display_error(SetupGUIController.main_window, f"Could not paste rows to this PypeIt file: {e}")


    def comment_out_metadata_rows(self):
        """Comment out one or more selected metadata rows."""
        if self._view is None:
            return
        row_indices = self._view.selectedRows()
        msgs.info(f"Commenting out {len(row_indices)} rows.")
        if len(row_indices) > 0:
            self._model.commentMetadataRows(row_indices)
    
    def uncomment_metadata_rows(self):
        """Uncomment previously commented out selected metadata rows."""
        if self._view is None:
            return
        row_indices = self._view.selectedRows()
        msgs.info(f"Unommenting out {len(row_indices)} rows.")
        if len(row_indices) > 0:
            self._model.uncommentMetadataRows(row_indices)

    def delete_metadata_rows(self) -> bool:
        """Remove one or more selected rows from the PypeItMetadataModel.
        
        Return: 
            True if there were metadata rows deleted, False if there
            weren't any rows to delete.
        """
        if self._view is None:
            return False
        row_indices = self._view.selectedRows()
        msgs.info(f"Removing {len(row_indices)} rows.")
        if len(row_indices) > 0:               
            self._model.removeMetadataRows(row_indices)
            return True
        return False


class PypeItObsLogController(QObject):
    """PypeItObsLogController responsible for responding to user interaction
    as part of a MVC triplet with PypeItObsLogView and PypeItObsLogModel

    Args:
        main_window (:obj:`UserPromptDelegate`): A view object that can prompt the user.
        model (:obj:`PypeItObsLogModel`): The model for the obs log.
        operation_thread (:obj:`pypeit.setup_gui.controller.SetupGUIController`): The main Setup GUI controller.
    """

    def __init__(self, model, setup_gui_controller):
        super().__init__()
        self._model = model
        self.setup_gui_controller = setup_gui_controller

    def setModel(self, model):
        self._model = model

    def getMetadataController(self, model):
        return PypeItMetadataController(model, is_pypeit_file=False)


    def setSpectrograph(self, spectrograph_name):
        self._model.set_spectrograph(spectrograph_name)

        if self._model.state != ModelState.NEW:
            # Re-run setup with the new spectrograph
            self.setup_gui_controller.run_setup()

        else:
            self._model.set_spectrograph(spectrograph_name)

    def removePaths(self, rows):
        # Remove paths in reverse order, so that indexes don't change when a row is removed
        for row in sorted(rows, reverse=True):
            self._model.paths_model.removeRow(row)

    def addNewPath(self, new_path):
        """Add a new path to the observation log"""
        msgs.info(f"Adding new path {new_path}")
        self._model.add_raw_data_directory(new_path)

class PypeItFileController(QObject):
    """PypeItFileController responsible for responding to user interaction
    as part of a MVC triplet with PypeItFileView and PypeItFileModel

    Args:
        main_window (:obj:`UserPromptDelegate`): A view object that can prompt the user.
        model (:obj:`PypeItFileModel`): The model for the obs log.
    """

    def __init__(self, model):
        self._model = model

    def setModel(self, model):
        self._model = model

    def getMetadataController(self, model):
        return PypeItMetadataController(model, is_pypeit_file=True)


class SetupGUIController(QObject):
    """Controller for the PypeIt setup gui. It is responsible for initializing the GUI,
    and performing actions requested by the user.
    
    Args:
        args (:class:`argparse.Namespace`): The non-QT command line arguments used to start the GUI.
    """


    main_window = None
    model = PypeItSetupGUIModel()

    def __init__(self, args):
        super().__init__()

        QCoreApplication.setOrganizationName("PypeIt")
        QCoreApplication.setApplicationName("SetupGUI")
        QCoreApplication.setOrganizationDomain("pypeit.readthedocs.io")

        if args.logfile is not None:
            logpath = Path(args.logfile)
            if logpath.exists():
                timestamp = datetime.now(__UTC__).strftime("%Y%m%d-%H%M%S")
                old_log=logpath.parent / (logpath.stem + f".{timestamp}" + logpath.suffix)
                logpath.rename(old_log)
                
        self.model.setup_logging(args.logfile, verbosity=args.verbosity)
        if args.spectrograph is not None:
            self.model.obslog_model.set_spectrograph(args.spectrograph)
        if args.root is not None:
            for root_dir in args.root:
                self.model.obslog_model.add_raw_data_directory(root_dir)
        if args.spectrograph is not None and args.root is not None:
            self.run_setup_at_startup = True
        else:
            self.run_setup_at_startup = False

        self.model.obslog_model.default_extension = args.extension
        SetupGUIController.main_window = SetupGUIMainWindow(self.model, self)
        self.operation_thread = OperationThread()


    def getObsLogController(self, model):
        """Create the PypeItObsLogController as part of a MVC triplet.

        Args:
            model (:obj:`PypeItObsLogModel`): The model for the obs log.
            main_window (:obj:`SetupGUIMainWindow`): The mainwindow of the setup GUI.
        """
        return PypeItObsLogController(model, self)
    
    def getPypeItFileController(self, model):
        """Create the PypeItObsLogController as part of a MVC triplet.

        Args:
            model (:obj:`PypeItFileModel`): The model for the obs log.
        """
        return PypeItFileController(model)

    def start(self, app):
        """
        Starts the PypeItSetupGUi event loop. Exits the GUI when the GUI is closed.

        Args:
            app (QApplication): The Qt application object for the GUI. The caller is expected
                                to pass any Qt specific command line arguments to this object
                                before calling start(). 
        """
        self.app = app
        self.main_window.show()
        if self.run_setup_at_startup:
            self.run_setup()

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
        sys.exit(app.exec_())

    def save_all(self):
        """"
        Save all unique configurations as pypeit files. Called in response to the user
        clicking the "Save All" button.
        """
        use_for_all_location = None
        for file_model in self.model.pypeit_files.values():                   

            # If the user clicked "Use this location for everthing"
            # Use that location for files that don't already have a location
            if file_model.save_location is None and use_for_all_location is not None:
                file_model.save_location = use_for_all_location

            # Save the tab, if neccessary prompting the user for a location to save to
            response = self._save_file(file_model, prompt_for_all=True)

            if response == DialogResponses.CANCEL:
                # Cancel was pressed
                return
            if response == DialogResponses.ACCEPT_FOR_ALL:
                # Use the selected location for everything going forward
                use_for_all_location = file_model.save_location


    def save_one(self):
        """ Saves the currently selected configuration as a pypeit file. Called in response to the user
        clicking the "Save Tab" button."""

        view_widget = self.main_window.tab_widget.currentWidget()
        if isinstance(view_widget, PypeItFileView):
            file_model = view_widget.model
            self._save_file(file_model)
        else:
            # Shouldn't really happen, it would mean the save tab button was enabled
            # when it shouldn't be. We'll handle this case and log it to prevent a crash
            # just in case though.
            msgs.warn(f"Attempt to save a tab that is *not* a PypeItFileView!")

    
    def _save_file(self, file_model : PypeItFileModel, prompt_for_all : bool=False) -> DialogResponses:
        """Helper method to save a file prompting the user for a location to save to
        if needed.
        
        Args:
            file_model:     The file model to save.
            prompt_for_all: Whether to prompt the user if they want to sue the location for
                            subsequent files. e.g. for a save_all operation.
        Return: 
            The DialogResponse from the user, or DialogResponses.ACCEPT if it wasn't
            neccessary to prompt the user.
        """
        msgs.info(f"Saving config {file_model.name_stem}")
        if file_model.save_location is None:
            dialog = FileDialog.create_save_location_dialog(self.main_window, file_model.name_stem, prompt_for_all=prompt_for_all)
            response = dialog.show()
            if response == DialogResponses.CANCEL:
                return response
            file_model.save_location = dialog.selected_path
        else:
            response = DialogResponses.ACCEPT

        try:
            file_model.save()
        except Exception as e:
            display_error(self.main_window, str(e))
        
        return response


    def clear(self):
        """Resets the GUI to it's initial state. Called in response to the user
        clicking the "Clear" button. This will prompt the user if they want to save
        any unsaved changes first."""

        # Prompt to save unsaved changes
        if self.model.state == ModelState.CHANGED:
            response = prompt_to_save(self.main_window)
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        self.model.reset()

    def close(self, file_model):

        if file_model.state == ModelState.CHANGED:
            # If the file has unsaved data ask the user if they want to save it
            response = prompt_to_save(self.main_window)
            if response == DialogResponses.SAVE:
                # Save the file, prompting for a location to save to if necessary
                self._save_file(file_model)
            elif response == DialogResponses.CANCEL:
                return False

        self.model.removeFile(file_model.name_stem)
        return True


    def exit(self):
        """Exits the GUI. Called in response to the user clicking the "Exit" button.
        This will prompt the user if they want to save any unsaved changes first."""

        # Prompt to save unsaved changes
        if self.model.state == ModelState.CHANGED:
            response = prompt_to_save(self.main_window)
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        self.app.quit()

    def run_setup(self):
        """Runs setup on the currently selected raw data directories. 
           Called in response to the user the clicking the "Run Setup" or changing the spectrograph.
           This will prompt the user if they want to save any unsaved changes first and
           start the operation in a background thread.
        """

        if self.model.state == ModelState.CHANGED:
            response = prompt_to_save(self.main_window)
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        self.operation_thread.startOperation(SetupOperation(self.model))

    def createNewPypeItFile(self):
        # First figure out the name
        # TODO, use New File convention as above, which forces "save as" when saving a new file,
        # or current default single letter naming scheme?
        # If we use letters, use base 26 numbers to support more than 1 letter?
        #new_name = "New File"
        #i = 1
        #while new_name in self.pypeit_files.keys():
        #    i += 1
        #    new_name = f"New File{i}"


        if len(self.model.pypeit_files) == 0:
            # No configs, just add "A"
            new_name = "A"
        else:
            largest_name = max(self.model.pypeit_files.keys())
            if largest_name == 'z':
                raise ValueError("Failed to create new setup because there are too many loaded.")
            new_name = chr(ord(largest_name)+1)

        self.model.createEmptyPypeItFile(new_name)


    def open_pypeit_file(self):
        """Opens a PypeIt file. Called in response to the user the clicking the "Open" button. 
           This method prompts the user to discard or save any current changes,
           prompts the user for a pypeit to open, and opens it.
        """
        # Prompt to save current changes if needed
        if self.model.state == ModelState.CHANGED:
            response = prompt_to_save(self.main_window)
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        open_dialog = FileDialog.create_open_file_dialog(self.main_window, "Select PypeIt File", file_type=FileType("PypeIt input files",".pypeit"))
        result = open_dialog.show()
        if result != DialogResponses.CANCEL:
            self.operation_thread.startOperation(OpenFileOperation(self.model, open_dialog.selected_path))
