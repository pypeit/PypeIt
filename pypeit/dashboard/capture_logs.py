# there is a few ways for me to capture logs

# one of the ways is to have access to the logger that pypeit is using and set up
# a handler that captures that output and puts it in a variable

# another way is to use subprocess and have capture output


from PyQt6.QtCore import QThread, pyqtSignal
import subprocess
import sys

class PypeitWorker(QThread):
    line_received = pyqtSignal(str)
    finished = pyqtSignal(int)
    def __init__(self,command):
        super().__init__()
        self.command = command
        self._stop_requested = False

    def run(self):
        process = subprocess.Popen(
                self.command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1 # line buffered
                )

        stdout = process.stdout
        if stdout is not None:
            for line in stdout:
                if self._stop_requested:
                    process.terminate() # terminates the process
                    break
                if line:
                    self.line_received.emit(line.rstrip())

        process.wait()
        self.finished.emit(process.returncode)

    def stop(self):
        self._stop_requested = True

