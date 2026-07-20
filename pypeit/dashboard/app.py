"""
Application bootstrap for the PypeIt Dashboard.

This mirrors the ``setup_gui`` startup path (``controller.start_gui`` +
``SetupGUIController.start``): it creates the :obj:`QApplication`, applies
the shared window-icon/font conventions, installs the loud-failure
excepthook, builds the main window from the reduction metadata, and runs
the Qt event loop.
"""

import sys
import signal

from qtpy.QtWidgets import QApplication
from qtpy.QtGui import QIcon
from qtpy.QtCore import QCoreApplication, QTimer

from pypeit import log, dataPaths
from pypeit.dashboard import util
from pypeit.dashboard.model import DashboardModel
from pypeit.dashboard.view.main_window import DashboardMainWindow

# The window icon shipped with the package (pypeit/data/images).
_ICON_PATH = dataPaths.images.path / 'pypeit_logo.png'


def _make_application():
    """
    Create the :obj:`QApplication` with the dashboard's icon and font
    conventions (reusing the ``setup_gui`` approach).

    Returns:
        ``QApplication``: The application instance.
    """
    app = QApplication(sys.argv)
    QCoreApplication.setOrganizationName('PypeIt')
    QCoreApplication.setApplicationName('Dashboard')
    QCoreApplication.setOrganizationDomain('pypeit.readthedocs.io')

    # Window/application icon (non-fatal if missing).
    if _ICON_PATH.is_file():
        app.setWindowIcon(QIcon(str(_ICON_PATH)))

    # Match setup_gui: ensure a comfortable minimum font size.
    font = app.font()
    if font.pointSizeF() < 14.0:
        font.setPointSize(14)
        app.setFont(font)

    return app


def launch(args):
    """
    Launch the PypeIt Dashboard and run the Qt event loop.

    Args:
        args (:obj:`argparse.Namespace`):
            Parsed command-line arguments.  Must provide ``pypeit_file``
            (the required positional ``.pypeit`` argument) and may provide
            ``redux_path``.

    Returns:
        int: The Qt application exit code.
    """
    # Make failures loud from the very first interaction (Debugging plan).
    util.install_excepthook()

    # Build the reduction-state model: load the state file if present, else
    # derive it (R4/R5; deriving may briefly block the UI on launch, which is
    # acceptable per the design).  Edge states (missing/empty/malformed) are
    # reported via model.load_status and rendered by the Status view (R11).
    redux_path = getattr(args, 'redux_path', None)
    model = DashboardModel(args.pypeit_file, redux_path=redux_path,
                           derive=True)
    if model.header_info is not None:
        log.info(f'Launching PypeIt Dashboard for '
                 f'{model.header_info.pypeit_file} '
                 f'({model.header_info.spectrograph}, '
                 f'{model.header_info.path}); state: {model.load_status}')
    else:
        log.warning(f'Launching PypeIt Dashboard: {model.load_status} '
                    f'({model.error})')

    app = _make_application()
    window = DashboardMainWindow(model)

    # Keep the window within the available screen area: the design's default
    # size (1650x900) is wider than some displays, so clamp to ~95% of the
    # available geometry and anchor at its top-left.
    screen = app.primaryScreen()
    if screen is not None:
        avail = screen.availableGeometry()
        width = min(window.width(), int(avail.width() * 0.95))
        height = min(window.height(), int(avail.height() * 0.95))
        window.resize(width, height)
        window.move(avail.left(), avail.top())

    window.show()

    # Allow Ctrl+C to interrupt the GUI from the terminal: Qt runs its event
    # loop in C, so we poke the Python interpreter on a timer (same trick as
    # setup_gui) to give the SIGINT handler a chance to run.
    def _sigint_handler(signalnum, frame):
        """
        Exit cleanly when Ctrl+C is pressed (the standard ``signal`` handler
        signature).

        Args:
            signalnum (:obj:`int`): The signal number received.
            frame: The interrupted stack frame (unused).
        """
        if signalnum == signal.SIGINT:
            log.info('Ctrl+C pressed; closing the PypeIt Dashboard.')
            sys.exit()

    signal.signal(signal.SIGINT, _sigint_handler)
    timer = QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)

    # qtpy exposes the Qt5-style exec_ as well as exec; use exec_ to match
    # setup_gui.
    return app.exec_()
