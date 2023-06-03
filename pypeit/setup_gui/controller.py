"""
The controller portion of the PypeIt Setup GUI.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import traceback
import signal
import sys
import datetime
import re
from pathlib import Path
from functools import wraps
from contextlib import contextmanager

from qtpy.QtCore import QCoreApplication, Signal
from qtpy.QtCore import QObject, Qt, QThread
from qtpy.QtGui import QKeySequence
from qtpy.QtWidgets import QAction

from pypeit.setup_gui.view import SetupGUIMainWindow, DialogResponses
from pypeit.setup_gui.model import PypeItSetupGUIModel, ModelState
from pypeit import msgs


class OpCanceledError(Exception):
    """Exception thrown when a background operation has been canceled."""
    def __init__(self):
        super().__init__()

class OperationThread(QThread):
    """Thread to run a background operation.
    
    Args:
        operation (:obj:`collections.abc.Callable`): The callable function or object that performs the operation.
    """

    completed = Signal(bool, tuple)

    def __init__(self, main_window):
        super().__init__()
        self._operation = None
        self._main_window = main_window
        self._max_progress = None

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

        if self._max_progress is None:
            if max_progress is not None:
                self._max_progress = max_progress
                self._main_window.create_progress_dialog(self._operation.name, max_progress, self._cancel_op)
            else:
                # Ignore the progress if there's no max progress yet
                return
            
        self._main_window.show_operation_progress(increase=1, 
                                                  message = progress_message)

    def _op_complete(self, canceled, exc_info):
        """Signal handler that is notified when a background operation completes.
        
        Args:
            canceled (bool): Whether or not the operation was canceled.
            exc_info (tuple): The exception information if the operation failed. None if it succeeded
        """
        msgs.info("Op complete")
        if self._operation is not None:
            self._operation.progressMade.disconnect(self._op_progress)
            self._main_window.operation_complete()
            self._operation.postRun(canceled, exc_info)
            self._operation = None
            self._max_progress = None
    

    def startOperation(self, operation):
        self._operation = operation
        if operation.preRun():
            operation.progressMade.connect(self._op_progress, type=Qt.QueuedConnection)
            self.completed.connect(self._op_complete, type=Qt.QueuedConnection)
            self.start()

class MetadataOperation(QObject):

    progressMade = Signal(int, str)

    def __init__(self, model, main_window):
        super().__init__()
        self._model=model
        self.main_window = main_window
        self.name = "Run Setup"
        self._max_progress = None

    def preRun(self):
        building_metadata_re = re.compile("Building metadata for (\d+) ")
        self._model.log_buffer.watch("building_metadata", building_metadata_re, self._buildingMetadata)
        
        added_metadata_re = re.compile("Adding metadata for (.*)$")
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
        self._model.log_buffer.unwatch("added_metadata")
        self._model.log_buffer.unwatch("building_metadata")

        if exc_info[0] is not None:
            traceback_string = "".join(traceback.format_exception(*exc_info))
            msgs.warn(f"Failed to {self.name.lower()}:\n" + traceback_string)
            self.main_window.display_error(f"Failed to {self.name.lower()} {exc_info[0]}: {exc_info[1]}")
            self._model.reset()
        elif canceled:
            self._model.reset()

    def run(self):
        pass

class SetupOperation(MetadataOperation):

    def __init__(self, model, main_window):
        super().__init__(model, main_window)

    def run(self):
        self._model.run_setup()

class OpenFileOperation(MetadataOperation):

    def __init__(self, model, main_window, file):
        super().__init__(model, main_window)
        self._file = file

    def run(self):
        self._model.open_pypeit_file(self._file)

class PypeItMetadataController(QObject):

    def __init__(self, main_window, model, is_pypeit_file):
        self._main_window = main_window
        self._model = model
        self.action_info = [("View File",   True,            None,                            self.view_file),
                            ("View Header", True,            None,                            self.view_header),
                            ("Copy",        True,            QKeySequence.StandardKey.Copy,   self.copy_metadata_rows),
                            ("Cut",         is_pypeit_file,  QKeySequence.StandardKey.Cut,    self.cut_metadata_rows),
                            ("Paste",       is_pypeit_file,  QKeySequence.StandardKey.Paste,  self.paste_metadata_rows),
                            ("Comment Out", is_pypeit_file,  None,                            self.comment_out_metadata_rows),
                            ("Uncomment",   is_pypeit_file,  None,                            self.uncomment_metadata_rows),
                            ("Delete",      is_pypeit_file,  QKeySequence.StandardKey.Delete, self.delete_metadata_rows) ]

    def getActions(self, parent):
        action_list = []
        for info in self.action_info:
            action = QAction(text=info[0], parent=parent)
            action.setEnabled(info[1])
            if info[2] is not None:
                action.setShortcut(info[2])
            action.triggered.connect(info[3])
            action_list.append(action)
        return action_list
            

    def view_file(self):
        pass

    def view_header(self):
        pass

    def copy_metadata_rows(self):
        pass

    def cut_metadata_rows(self):
        pass

    def paste_metadata_rows(self):
        pass

    def comment_out_metadata_rows(self):
        pass
    
    def uncomment_metadata_rows(self):
        pass

    def delete_metadata_rows(self):
        pass



class PypeItObsLogController(QObject):
    """PypeItObsLogController responsible for responding to user interaction
    as part of a MVC triplet with PypeItObsLogView and PypeItObsLogModel

    Args:
        main_window (:obj:`UserPromptDelegate`): A view object that can prompt the user.
        model (:obj:`PypeItObsLogModel`): The model for the obs log.
        operation_thread (:obj:`pypeit.setup_gui.controller.SetupGUIController`): The main Setup GUI controller.
    """

    def __init__(self, main_window, model, setup_gui_controller):
        super().__init__()
        self._model = model
        self._main_window = main_window
        self.setup_gui_controller = setup_gui_controller

    def setModel(self, model):
        self._model = model

    def getMetadataController(self, model):
        return PypeItMetadataController(self._main_window, model, is_pypeit_file=False)


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

    def __init__(self, main_window, model):
        self._model = model
        self._main_window = main_window

    def setModel(self, model):
        self._model = model

    def getMetadataController(self, model):
        return PypeItMetadataController(self._main_window, model, is_pypeit_file=True)


class SetupGUIController(QObject):
    """Controller for the PypeIt setup gui. It is responsible for initializing the GUI,
    and performing actions requested by the user.
    
    Args:
        args (:class:`argparse.Namespace`): The non-QT command line arguments used to start the GUI.
    """

    def __init__(self, args):
        super().__init__()

        QCoreApplication.setOrganizationName("PypeIt")
        QCoreApplication.setApplicationName("SetupGUI")
        QCoreApplication.setOrganizationDomain("pypeit.readthedocs.io")

        if args.logfile is not None:
            logpath = Path(args.logfile)
            if logpath.exists():
                timestamp = datetime.datetime.utcnow().strftime("%Y%m%d-%H%M%S")
                old_log=logpath.parent / (logpath.stem + f".{timestamp}" + logpath.suffix)
                logpath.rename(old_log)
                
        self.model = PypeItSetupGUIModel()
        self.model.setup_logging(args.logfile, verbosity=args.verbosity)
        if args.spectrograph is not None:
            self.model.obslog_model.set_spectrograph(args.spectrograph)
        if args.root is not None:
            self.model.obslog_model.add_raw_data_directory(args.root)
        if args.spectrograph is not None and args.root is not None:
            self.run_setup_at_startup = True
        else:
            self.run_setup_at_startup = False

        self.model.obslog_model.default_extension = args.extension
        self.main_window = SetupGUIMainWindow(self.model, self)
        self.operation_thread = OperationThread(self.main_window)

    def getObsLogController(self, model, main_window):
        """Create the PypeItObsLogController as part of a MVC triplet.

        Args:
            model (:obj:`PypeItObsLogModel`): The model for the obs log.
            main_window (:obj:`SetupGUIMainWindow`): The mainwindow of the setup GUI.
        """
        return PypeItObsLogController(main_window, model, self)
    
    def getPypeItFileController(self, model):
        """Create the PypeItObsLogController as part of a MVC triplet.

        Args:
            model (:obj:`PypeItFileModel`): The model for the obs log.
        """
        return PypeItFileController(self.main_window, model)

    def start(self, app):
        """
        Starts the PypeItSetupGUi event loop. Exits the GUI when the GUI is closed.

        Args:
            app (QApplication): The Qt application object for the GUI. The caller is expected
                                to pass any Qt specific command line arguments to this object
                                before calling start(). 
        """
        self.main_window.show()
        if self.run_setup_at_startup:
            self.run_setup()

        # QT Gobbles up the Python Ctrl+C handling, so the default PypeIt Ctrl+C handler won't
        # be called. So we reset it to the OS default
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        sys.exit(app.exec_())

    def save_all(self):
        """"
        Save all unique configurations as pypeit files. Called in response to the user
        clicking the "Save All" button.
        """
        try:
            response = DialogResponses.CANCEL
            location = None
            for file_model in self.model.pypeit_files.values():
                if file_model.save_location is None:
                    if location is None or response != DialogResponses.ACCEPT_FOR_ALL:
                        (location, response) = self.main_window.prompt_for_save_location(file_model.config_name, True)
                        if response == DialogResponses.CANCEL:
                            # Cancel was pressed
                            return

                file_model.save_location = location
                file_model.save()

        except Exception as e:
            self.view.display_error(str(e))

    def save_one(self):
        """ Saves the currently selected configuration as a pypeit file. Called in response to the user
        clicking the "Save Tab" button."""
        try:
            config_name = self.main_window.tab_widget.currentWidget().name
        except Exception as e:
            self.main_window.display_error(str(e))

        self._save_tab(config_name)
    
    def _save_tab(self, config_name):
        msgs.info(f"Saving config {config_name}")
        file_model = self.model.pypeit_files[config_name]
        if file_model.save_location is None:
            location, response = self.main_window.prompt_for_save_location(config_name)
            if response == DialogResponses.CANCEL:
                return
            else:
                file_model.save_location = location
        file_model.save()

    def clear(self):
        """Resets the GUI to it's initial state. Called in response to the user
        clicking the "Clear" button. This will prompt the user if they want to save
        any unsaved changes first."""

        # Prompt to save unsaved changes
        if self.model.state == ModelState.CHANGED:
            response = self.main_window.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        self.model.reset()

    def createNewFile(self, source_file_name, selectedRows):
        try:
            self.model.createNewPypeItFile(source_file_name, selectedRows)
        except Exception as e:
            self.main_window.display_error(f"Failed to create new tab {e.__class__.__name__}: {e}")
            msgs.warn(f"Failed to create new tab.")
            msgs.warn(traceback.format_exc())

    def closeFile(self, file_name):
        file_model = self.model.pypeit_files[file_name]
        if file_model.state == ModelState.CHANGED:
            response = self.main_window.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_one(file_name)
            elif response == DialogResponses.CANCEL:
                return
        self.model.removeConfig(file_name)

    def exit(self):
        """Exits the GUI. Called in response to the user clicking the "Exit" button.
        This will prompt the user if they want to save any unsaved changes first."""

        # Prompt to save unsaved changes
        if self.model.state == ModelState.CHANGED:
            response = self.main_window.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        sys.exit(0)

    def run_setup(self):
        """Runs setup on the currently selected raw data directories. 
           Called in response to the user the clicking the "Run Setup" or changing the spectrograph.
           This will prompt the user if they want to save any unsaved changes first and
           start the operation in a background thread.
        """

        if self.model.state == ModelState.CHANGED:
            response = self.main_window.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        self.operation_thread.startOperation(SetupOperation(self.model, self.main_window))

    def open_pypeit_file(self):
        """Opens a PypeIt file. Called in response to the user the clicking the "Open" button. 
           This method prompts the user to discard or save any current changes,
           prompts the user for a pypeit to open, and opens it.
        """
        # Prompt to save current changes if needed
        if self.model.state == ModelState.CHANGED:
            response = self.main_window.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        file_to_open = self.main_window.prompt_for_open_file("Select PypeIt File", filter="PypeIt input files (*.pypeit)")

        if file_to_open is not None:
            self.operation_thread.startOperation(OpenFileOperation(self.model, self.main_window, file_to_open))
