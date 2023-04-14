"""
The controller portion of the PypeIt Setup GUI.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import traceback
import signal
import sys
import datetime
from pathlib import Path

from qtpy.QtCore import QCoreApplication
from qtpy.QtCore import QObject, Qt, QThread


from pypeit.setup_gui.view import SetupGUIMainWindow, DialogResponses
from pypeit.setup_gui.model import PypeItSetupModel, ModelState
from pypeit import msgs


class OperationThread(QThread):
    """Thread to run a background operation.
    
    Args:
        operation (:obj:`collections.abc.Callable`): The callable function or object that performs the operation."""

    def __init__(self, operation):
        super().__init__()
        self._operation = operation

    def run(self):
        """Runs an operation in a background thread."""
        self._operation()


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
                
        self.model = PypeItSetupModel()
        self.model.setup_logging(args.logfile, verbosity=args.verbosity)
        if args.spectrograph is not None:
            self.model.set_spectrograph(args.spectrograph)
        if args.root is not None:
            self.model.add_raw_data_directory(args.root)
        if args.spectrograph is not None and args.root is not None:
            self.run_setup_at_startup = True
        else:
            self.run_setup_at_startup = False

        self.model.default_extension = args.extension
        self.view = SetupGUIMainWindow(self.model, self)
        self._thread = None


    def start(self, app):
        """
        Starts the PypeItSetupGUi event loop. Exits the GUI when the GUI is closed.

        Args:
            app (QApplication): The Qt application object for the GUI. The caller is expected
                                to pass any Qt specific command line arguments to this object
                                before calling start(). 
        """
        self.view.show()
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
                        (location, response) = self.view.prompt_for_save_location(file_model.config_name, True)
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
            config_name = self.view.setup_view.currentWidget().name
        except Exception as e:
            self.view.display_error(str(e))

        self._save_tab(config_name)
    
    def _save_tab(self, config_name):
        msgs.info(f"Saving config {config_name}")
        file_model = self.model.pypeit_files[config_name]
        if file_model.save_location is None:
            location, response = self.view.prompt_for_save_location(config_name)
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
            response = self.view.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        self.model.reset()

    def removeConfig(self, config_name):
        file_model = self.model.pypeit_files[config_name]
        if file_model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_one(config_name)
            elif response == DialogResponses.CANCEL:
                return
        self.model.removeConfig(config_name)

    def exit(self):
        """Exits the GUI. Called in response to the user clicking the "Exit" button.
        This will prompt the user if they want to save any unsaved changes first."""

        # Prompt to save unsaved changes
        if self.model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        sys.exit(0)

    def open_pypeit_file(self):
        """Opens a PypeIt file. Called in response to the user the clicking the "Open" button. 
           This method prompts the user to discard or save any current changes,
           prompts the user for a pypeit to open, and opens it.
        """
        # Prompt to save current changes if needed
        if self.model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        file_to_open = self.view.prompt_for_open_file("Select PypeIt File", filter="PypeIt input files (*.pypeit)")

        if file_to_open is not None:
            try:
                self.model.open_pypeit_file(file_to_open)
            except Exception as e:
                msgs.warn(f"Failed to open {file_to_open} : {e}")
                msgs.warn(traceback.format_exc())
                self.view.display_error(f"Failed to open {file_to_open}\n{e}")
                self.model.reset()

    def run_setup(self):
        """Runs setup on the currently selected raw data directories. 
           Called in response to the user the clicking the "Run Setup" button. 
           This will prompt the user if they want to save any unsaved changes first and
           start the operation in a background thread.
        """

        if self.model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        num_raw_data_files =  self.model.scan_raw_data_directories()
        if num_raw_data_files == 0:
            self.view.display_error(self.tr("Could not find any files in any of the raw data paths."))
            return

        # Start the setup operation in a different thread
        self.start_operation("Reading Files...", self.model.run_setup, self.run_setup_complete, num_raw_data_files)

    def start_operation(self, op_name, op_function, op_cleanup_func, max_progress_value):
        """Start a background operation.  The operation runs in the background and displays a progress
        dialog to the user.
        
        Args:
            op_name (str):            The name of the operation to display the user.
            max_progress_value (int): The maximum progress value for the progress dialog. The progress
                                      always starts at 0 and is increased until it reaches this maximum value.
            op_function (:class:`collections.abc.Callable`):     The callable function or object to perform the operation.
            op_cleanup_func (:class:`collections.abc.Callable`): The callable function or object responsible for cleaning up after
                                                                 the operation.
        """
        
        self.view.start_operation_progress(op_name, max_progress_value, self._cancel_op)

        # Attach to model's signals, using a queued connectio so that the model can run in its own thread
        self.model.operation_progress.connect(self._op_progress, type=Qt.QueuedConnection)
        self.model.operation_complete.connect(op_cleanup_func, type=Qt.QueuedConnection)

        # Start the thread to generate an obslog from raw data
        self._thread = OperationThread(op_function)
        self._thread.start()

    def _cancel_op(self):
        """Cancels an in progress background operation when the user cancels the progress dialog."""
        if self._thread is not None:
            self._thread.requestInterruption()

    def _op_progress(self, progress_message=None):
        """Signal handler that receives progress information from the model as a background
        operation proceeds. It passes this to the view to increase the value showing in the
        progress dialog."""
        msgs.info(f"New Progress: {progress_message}")

        self.view.increase_operation_progress(increase=1, message = progress_message)

    def _op_complete(self, error_message):
        """Signal handler that is notified when a background operation completes.
        
        Args:
            error_message (str): An error message if the operation failed, "CANCEL" if the
                                 operation was canceled, or a blank string if the operation 
                                 succeeded.
        """
        msgs.info("Op complete")
        self._thread.wait()
        self.model.operation_progress.disconnect(self._op_progress)
        self.view.operation_complete()
        self._thread = None
        if error_message is not None and len(error_message) > 0:
            if error_message != "CANCEL":
                self.view.display_error(error_message)
            self.model.reset()

    def run_setup_complete(self, error_message):
        """Cleans up after setup has been run.

        Args:
            error_message (str): An error message if the operation failed, "CANCEL" if the 
                                 operation was canceled, or a blank string if the operation 
                                 succeeded.
        """
        self._op_complete(error_message)
        self.model.operation_complete.disconnect(self.run_setup_complete)