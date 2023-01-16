import traceback
import signal
import sys
import datetime

from qtpy.QtCore import QCoreApplication
from qtpy.QtWidgets import QApplication
from qtpy.QtCore import QObject, Qt, QThread, QSettings


from pypeit.pypmsgs import PypeItError
from pypeit.setup_gui.view import SetupGUIMainWindow, DialogResponses
from pypeit.setup_gui.model import PypeItSetupProxy, ModelState
from pypeit import msgs
from pathlib import Path

class OperationThread(QThread):

    def __init__(self, operation):
        super().__init__()
        self._operation = operation

    def run(self):
        self._operation()


class SetupGUIController(QObject):


    def __init__(self, args):
        super().__init__()

        QCoreApplication.setOrganizationName("PypeIt")
        QCoreApplication.setApplicationName("SetupGUI")
        QCoreApplication.setOrganizationDomain("pypeit.readthedocs.io")

        timestamp = datetime.datetime.utcnow().strftime("%Y%m%d-%H%M")
        logname = f"setup_gui_{timestamp}.log"

        self.model = PypeItSetupProxy()
        self.model.setup_logging(logname, args.verbosity)
        self.view = SetupGUIMainWindow(self.model, self)
        self._thread = None


    def start(self, app):
        """
        Start the PypeItSetupGUi
        """
        self.view.show()

        # QT Gobbles up the Python Ctrl+C handling, so the default PypeIt Ctrl+C handler won't
        # be called. So we reset it to the OS default
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        sys.exit(app.exec_())

    def save_all(self, location):
        try:
            self.model.save_all(location)
        except Exception as e:
            self.view.display_error(str(e))


    def save_tab(self, tab_name, location):    
        try:
            self.model.save_config(tab_name, location)
        except Exception as e:
            self.view.display_error(str(e))

    def clear(self):

        # Prompt to save unsaved changes
        if self.model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        self.model.reset()

    def exit(self):

        # Prompt to save unsaved changes
        if self.model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        sys.exit(0)

    def open_pypeit_file(self):

        # Prompt to save current changes if needed
        if self.model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        # Get history from settings
        settings = QSettings()
        settings.beginGroup("OpenPypeItFile")
        history = settings.value("History")
        if history is None:
            history = []
        elif isinstance(history, str):
            history = [history]

        file_to_open = self.view.prompt_for_file("Select PypeIt File", filter="PypeIt input files (*.pypeit)", history = history)
        if file_to_open is not None:
            # The FileDialog only wants paths in its history
            path = str(Path(file_to_open).parent)
            if path not in history:
                history.append(path)
                settings.setValue("History", history)

            try:
                self.model.open_pypeit_file(file_to_open)
            except Exception as e:
                msgs.warn(f"Failed to open {file_to_open} : {e}")
                msgs.warn(traceback.format_exc())
                self.view.display_error(f"Failed to open {file_to_open}\n{e}")
                self.model.reset()


    def start_operation(self, op_name, op_function, op_cleanup_func, max_progress_value):
        
        self.view.start_operation_progress(op_name, max_progress_value, self._cancel_op)

        # Attach to model's signals, using a queued connectio so that the model can run in its own thread
        self.model.operation_progress.connect(self._op_progress, type=Qt.QueuedConnection)
        self.model.operation_complete.connect(op_cleanup_func, type=Qt.QueuedConnection)

        # Start the thread to generate an obslog from raw data
        self._thread = OperationThread(op_function)
        self._thread.start()

    def _cancel_op(self):
        if self._thread is not None:
            self._thread.requestInterruption()

    def _op_progress(self, progress_message=None):
        msgs.info(f"New Progress: {progress_message}")

        self.view.increase_operation_progress(increase=1, message = progress_message)

    def _op_complete(self, error_message):
        msgs.info("Op complete")
        self._thread.wait()
        self.model.operation_progress.disconnect(self._op_progress)
        self.view.operation_complete()
        self._thread = None
        if error_message is not None and len(error_message) > 0:
            if error_message != "CANCEL":
                self.view.display_error(error_message)
            self.model.reset()

    def spectrograph_changed(self, new_spec):
        if new_spec == self.model.spectrograph:
            return

        msgs.info(f"New spectrograph '{repr(new_spec)}' new index: {self.view.spectrograph.currentIndex()}")

        if self.model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        self.model.reset()
        if new_spec in self.model.available_spectrographs:
            self.model.set_spectrograph(new_spec)
        else:
            self.model.set_spectrograph(None)

    def run_setup(self, new_raw_data_directory):
        if new_raw_data_directory is None or not isinstance(new_raw_data_directory, str) or new_raw_data_directory == "":
            # We're re-running with the old directory
            new_raw_data_directory = self.model.raw_data_directory

        if self.model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        try:
            self.model.set_raw_data_directory(new_raw_data_directory)
        except PypeItError:
            self.view.display_error(f"Could not find any files in {new_raw_data_directory}")
            return

        num_raw_data_files =  self.model.get_num_raw_files()
        
        # Start the setup operation in a different thread
        self.start_operation("Reading Files...", self.model.generate_obslog, self.generate_obslog_complete, num_raw_data_files)

    def generate_obslog_complete(self, error_message):
        self._op_complete(error_message)
        self.model.operation_complete.disconnect(self.generate_obslog_complete)