"""Classes for displaying text content to a Qt window.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pathlib import Path
import io
from typing import Optional,Union

from qtpy.QtWidgets import  QHBoxLayout, QVBoxLayout, QFileDialog, QWidget, QPlainTextEdit, QPushButton

from qtpy.QtGui import QIcon, QFont,QTextCursor
from qtpy.QtCore import Qt, Signal, QSettings, QEvent

from pypeit import msgs
from pypeit.setup_gui.dialog_helpers import display_error, FileDialog, FileType, PersistentStringListModel, DialogResponses
from pypeit.setup_gui.model import LogBuffer

class TextViewerWindow(QWidget):
    """Window to display text. The text can also be saved to a file.
    
    Args:
        title:        Title for the window
        width:        Initial width of the window, in characters.
        height:       Initial height of the window, in characters.
        filename:     Name of the file being viewed.
        text_stream:  The file object for the text to display.
        start_at_top: Whether to start with the window scrolled to the top of the text.
        filters:      File extension filters for saving the file, in the format used by QFileDialog.
                      Defaults to ["Text Files ('*'.txt)"]
    """

    closed = Signal()
    """Signal sent when the window is closed."""


    def __init__(self, title : str, width : int, height : int, text_stream : io.TextIOBase, start_at_top : bool, filename: Optional[Union[str,Path]] = None, file_type : FileType = FileType("Text Files", ".txt")):
        super().__init__()
        self._text_stream = text_stream
        self._file_type = file_type

        # TODO , to be a more general purpose text viewer this could open the file if no text_stream is given,
        # but I don't need that ability right now, this is just used for saving
        self._filename = filename

        self.setWindowIcon(QIcon(str(Path(__file__).parent / "images/window_icon.png")))
        self.setWindowTitle(title)
        logLayout = QVBoxLayout(self)

        # Get a fixed width font
        fixed_font = QFont()
        fixed_font.setFamilies(["Monospace", "Courier"])
        fixed_font.setFixedPitch(True)
        fixed_font.setPointSize(12)

        # Create the log viewer widget
        self.textViewer=QPlainTextEdit(parent=self)
        self.textViewer.setFont(fixed_font)
        self.textViewer.setReadOnly(True)

        # Set the window size for widthxheight characters
        viewer_margins = self.textViewer.contentsMargins()
        layout_margins = logLayout.contentsMargins()
        parent_margins = self.contentsMargins()
        font_metrics = self.textViewer.fontMetrics()
        char_width = font_metrics.averageCharWidth()*width
        char_height = font_metrics.height()*height
        
        self.resize(parent_margins.left() + parent_margins.right() +
                    layout_margins.left() + layout_margins.right() +
                    viewer_margins.left() + viewer_margins.right() +
                    char_width, 
                    parent_margins.top() + parent_margins.bottom() +
                    layout_margins.top() + layout_margins.bottom() +
                    viewer_margins.top() + viewer_margins.bottom() +
                    char_height)

        logLayout.addWidget(self.textViewer)

        # Fill the viewer with the contents from the text stream
        for line in self._text_stream:
            self.addText(line)

        # Scroll to the top after adding all of the text
        if start_at_top:
            self.textViewer.moveCursor(QTextCursor.MoveOperation.Start)

        # Buttons
        buttonLayout = QHBoxLayout()
        logLayout.addLayout(buttonLayout)

        # Button to save the log
        self.saveButton = QPushButton(text=self.tr("Save..."))
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
        history = PersistentStringListModel("TextDirectory", "History")

        # Create the dialog.
        save_dialog = FileDialog(self, self.tr(f"Enter file name"), QFileDialog.AnyFile, default_file=self._filename,
                                 history=history, file_type=self._file_type, save=True)

        # Show the dialog.
        response = save_dialog.show()
        if response  != DialogResponses.CANCEL:            
            # Save the file
            try:
                with open(save_dialog.selected_path, "w") as f:
                    if self._text_stream.seekable():
                        self._text_stream.seek(0)
                    for message in self._text_stream:
                        f.write(message)
                msgs.info(f"File saved to {save_dialog.selected_path}.")
                self._filename =  save_dialog.selected_path
            except Exception as e:
                display_error(self, str(e))

    def closeEvent(self, event):
        """Event handle for closing the log window. Overridden from QWidget.

        Args:
            event (QEvent): The close event. Not used by this implementation.
        """
        # Notify the main window.
        self.closed.emit()
        return super().closeEvent(event)

    def addText(self, text: str):
        """Add text to the text viewer
        
        Args:
            text: The text to add
        """

        # Add the message to our internal text document
        self.textViewer.appendPlainText(text.strip())

class LogWindow(TextViewerWindow):
    """ Text window for viewing the PypeIt log messages as they occurr.
        logBuffer:   LogBuffer receiving PypeIt log messages.
    """
    newMessage = Signal(str)
    """Signal sent when a new log message is logged."""

    def __init__(self, logBuffer : LogBuffer):
        super().__init__(text_stream=logBuffer,
                         width=100, height=50,
                         start_at_top=False,
                         title = "PypeIt Setup GUI Log",
                         file_type = FileType("Log Files", ".log"))
    
        self._logBuffer = logBuffer

        # Preallocate the maximum length of the log buffer
        self.textViewer.document().setMaximumBlockCount(self._logBuffer.maxlen)

        # Listen for new log messages and append them to the viewer
        # This is done with a callbackto self._messageLogged, which then
        # emits an event from self.newMessage to notify self._addMessage.
        # This is done so that log messages from a different thread are
        # properly queued to the GUIs event thread.
        self.newMessage.connect(self.addText,Qt.QueuedConnection)
        self._logBuffer.watch("log window", None, self._messageLogged)

    def closeEvent(self, event: QEvent):
        """Event handle for closing the log window. Overridden from QWidget.

        Args:
            event: The close event. Not used by this implementation.
        """
        # Stop watching the log, and notify the main window.
        self.newMessage.disconnect(self.addText)
        self._logBuffer.unwatch("log window")
        return super().closeEvent(event)

    def _messageLogged(self, message : str):
        """Callback that is notified by the log buffer for each message logged.

        Args:
            message: The log message
        """
        # Emit the message so it can be picked up in the GUI's event thread.
        self.newMessage.emit(message)
