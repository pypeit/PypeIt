import traceback

from qtpy.QtCore import QObject, Qt, QThread, QSettings
from qtpy.QtWidgets import QApplication

from pypeit.pypmsgs import PypeItError
from pypeit.setup_gui.view import PypeItSetupView
from pypeit.setup_gui.model import ModelState
from pypeit import msgs
from pathlib import Path

class OperationThread(QThread):

    def __init__(self, operation):
        super().__init__()
        self._operation = operation

    def run(self):
        self._operation()


class PypeItSetupController(QObject):


    def __init__(self, model, view):
        super().__init__()

        self.view = view
        self.model = model
        self._thread = None
        self.view.setModel(model)
        self._connect_signals()        
        self.update_state()


    def _connect_signals(self):
        self.view.raw_data.location.activated.connect(self.raw_data_changed)
        self.view.spectrograph.currentTextChanged.connect(self.spectrograph_changed)
        self.view.exitButton.clicked.connect(QApplication.exit)
        self.view.setupButton.clicked.connect(self.run_setup)
        self.view.saveAllButton.clicked.connect(self.save_all)
        self.view.saveTabButton.clicked.connect(self.save_tab)
        self.view.openButton.clicked.connect(self.open_pypeit_file)
        self.view.clearButton.clicked.connect(self.clear)

    def update_state(self):

        if self.model.state == ModelState.NEW:
            self.view.spectrograph.setEnabled(True)
            if self.view.spectrograph.currentIndex() == -1:
                self.view.raw_data.setEnabled(False)
            else:
                self.view.raw_data.setEnabled(True)
            self.view.clearButton.setEnabled(False)
            self.view.setupButton.setEnabled(False)
            self.view.saveTabButton.setEnabled(False)
            self.view.saveAllButton.setEnabled(False)
        elif self.model.state == ModelState.CHANGED:
            self.view.clearButton.setEnabled(True)
            self.view.setupButton.setEnabled(True)
            self.view.saveAllButton.setEnabled(True)
        else: # Unchanged
            self.view.clearButton.setEnabled(True)
            self.view.setupButton.setEnabled(True)
            self.view.spectrograph.setEnabled(True)
            self.view.raw_data.setEnabled(True)
            self.view.saveTabButton.setEnabled(False)
            self.view.saveAllButton.setEnabled(False)

    def save_all(self):
        location = self.view.outdir.current_location
        if location is None or location == "":
            location=self.view.outdir.browse()
            if location is None:
                # Prompt was canceled
                return
        try:
            self.model.save_all(location)
        except Exception as e:
            self.view.display_error(str(e))

        self.update_state()

    def save_tab(self):
        tab_to_save = self.view.get_current_tab_name()
        if tab_to_save == "ObsLog":
            msgs.bug("Save tab pressed on ObsLog? That shouldn't happen...")
            return
    
        location = self.view.outdir.current_location
        if location is None or location == "":
            location=self.view.outdir.browse()
            if location is None:
                # Prompt was canceled
                return

        self.model.save_config(tab_to_save, location)
        self.view.saveTabButton.setEnabled(False)

    def clear(self):

        # Prompt to save unsaved changes
        if self.model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == PypeItSetupView.Responses.SAVE:
                self.save_all()
            elif response == PypeItSetupView.Responses.CANCEL:
                return

        self.model.reset()


    def open_pypeit_file(self):

        # Prompt to save current changes if needed
        if self.model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == PypeItSetupView.Responses.SAVE:
                self.save_all()
            elif response == PypeItSetupView.Responses.CANCEL:
                return

        # Get history from settings
        settings = QSettings()
        settings.beginGroup("OpenPypeItFile")
        history = settings.value("History")
        if history is None:
            history = []
        elif isinstance(history, str):
            history = [history]

        file_to_open = self.view.promptForFile("Select PypeIt File", filter="PypeIt input files (*.pypeit)", history = history)
        if file_to_open is not None:
            # The FileDialog only wants paths in its history
            path = str(Path(file_to_open).parent)
            if path not in history:
                history.append(path)
                settings.setValue("History", history)

            if self.model.state == ModelState.CHANGED:
                response = self.view.prompt_for_save()
                if response == PypeItSetupView.Responses.SAVE:
                    self.save_all()
                elif response == PypeItSetupView.Responses.CANCEL:
                    return

            try:
                self.model.open_pypeit_file(file_to_open)
                self.view.spectrograph.setCurrentText(self.model.spectrograph)
                self.view.createConfigTabs(self.model.configs.values())
            except Exception as e:
                msgs.warn(f"Failed to open {file_to_open} : {e}")
                msgs.warn(traceback.format_exc())
                self.view.display_error(f"Failed to open {file_to_open}\n{e}")
                self.model.reset()
            self.update_state()


    def start_operation(self, op_name, op_function, op_cleanup_func, max_progress_value):
        
        self.view.start_operation_progress(op_name, max_progress_value, self._cancel_op)

        # Attach to model's signals, using a queued connectio so that the model can run in its own thread
        self.model.operation_progress.connect(self._op_progress, type=Qt.QueuedConnection)
        self.model.operation_complete.connect(op_cleanup_func, type=Qt.QueuedConnection)

        # Start the thread to generate an obslog from raw data
        self._thread = OperationThread(op_function)
        self._thread.start()




    def raw_data_changed(self, new_location):
        msgs.info(f"New location '{repr(new_location)}' new text: {self.view.raw_data.location.currentText()}")
        
        new_raw_data_directory = self.view.raw_data.location.currentText()
        if self.model.raw_data_directory != new_raw_data_directory:
            self.run_setup(new_raw_data_directory)


    def _cancel_op(self):
        if self._thread is not None:
            self._thread.requestInterruption()

    def _op_progress(self, progress_message=None):
        msgs.info(f"New Progress: {progress_message}")

        self.view.set_operation_progress(increase=1, message = progress_message)

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
        self.update_state()

    def spectrograph_changed(self, new_spec):
        if new_spec == self.model.spectrograph:
            return

        msgs.info(f"New spectrograph '{repr(new_spec)}' new index: {self.view.spectrograph.currentIndex()}")

        if self.model.state == ModelState.UNCHANGED:
            self.model.reset()
        elif self.model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == PypeItSetupView.Responses.SAVE:
                self.save_all()
            elif response == PypeItSetupView.Responses.CANCEL:
                return
            self.model.reset()

        if new_spec in self.model.available_spectrographs:
            self.model.set_spectrograph(new_spec)
            self.view.raw_data.setEnabled(True)
            self.view.outdir.setEnabled(True)
        else:
            self.model.set_spectrograph(None)
            self.view.spectrograph.setEnabled(True)
            self.view.raw_data.setEnabled(False)
            self.view.outdir.setEnabled(False)

    def run_setup(self, new_raw_data_directory=None):
        if new_raw_data_directory is None or not isinstance(new_raw_data_directory, str) or new_raw_data_directory == "":
            # We're re-running with the old directory
            new_raw_data_directory = self.model.raw_data_directory

        if self.model.state == ModelState.CHANGED:
            response = self.view.prompt_for_save()
            if response == PypeItSetupView.Responses.SAVE:
                self.save_all()
            elif response == PypeItSetupView.Responses.CANCEL:
                return

        self.model.reset()
        try:
            self.model.set_raw_data_directory(new_raw_data_directory)
        except PypeItError:
            self.view.display_error(f"Could not find any files in {new_raw_data_directory}")
            return

        num_raw_data_files =  self.model.get_num_raw_files()
        
        # Get progress dialog from view
        self.start_operation("Reading Files...", self.model.generate_obslog, self.generate_obslog_complete, num_raw_data_files)

    def generate_obslog_complete(self, error_message):
        self._op_complete(error_message)
        self.model.operation_complete.disconnect(self.generate_obslog_complete)
        if error_message is None or len(error_message) == 0:
            self.view.createConfigTabs(self.model.configs.values())
