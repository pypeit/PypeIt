"""
PypeIt logging

Implementation heavily references loggers from astropy and sdsstools.
"""

import copy
import inspect
import io
import logging
from pathlib import Path
import re
import sys
import traceback
import warnings

from IPython import embed

# NOTE: BEWARE of importing anything from pypeit into this module.  It is likely
# to cause a circular import.

# Add a logging TEST level between INFO and WARNING
TEST_LEVEL_NUM = 25
logging.addLevelName(TEST_LEVEL_NUM, "TEST")
def logtest(self, message, *args, **kwargs):
    """Define TEST log level for calling"""
    # Skip this wrapper frame so filename/lineno point at the caller of log.test()
    kwargs.setdefault("stacklevel", 2)
    self.log(TEST_LEVEL_NUM, message, *args, **kwargs)
logging.Logger.test = logtest

# TODO: Can we put this *inside* the logger?
def short_warning(message, category, filename, lineno, line=None):
    """
    Overrides default formatting of warning messages.  The only arguments used
    are ``message`` and ``category``.  See :func:`warnings.formatwarning`.
    """
    return f'{category.__name__}: {message}'
warnings.formatwarning = short_warning

# NOTE: This is essentially a hack to deal with all the RankWarnings that numpy
# can throw during polynomial fitting.  Specifically this happens frequently
# in pypeit.core.fitting.PypeItFit.fit.  We should instead determine why these
# rank warnings are happening and address the root cause!
# 'default' means: "print the first occurrence of matching warnings for each
# location (module + line number) where the warning is issued"
# See: https://docs.python.org/3/library/warnings.html#warning-filter
import numpy as np
warnings.simplefilter('default', np.exceptions.RankWarning)


def color_text(
    text:str,
    color:list[int],
    bold:bool = False,
    nchar:int | None = None
) -> str:
    """
    Return an input string with escape characters to colorize text written to
    consoles.

    Parameters
    ----------
    text
        Text to colorize
    color
        3-element list of integers with the RGB color values
    bold
        Flag to make the text bold
    nchar
        Force the output text to be right-justified with this number of
        characters

    Returns
    -------
        Reformatted string
    """
    msg = '\033[1;' if bold else '\033['
    _text = f'{text}' if nchar is None else f'{text:>{nchar}}'
    return f'{msg}38;2;{color[0]};{color[1]};{color[2]}m{_text}\033[0m'


def clear_text_color(text:str) -> str:
    """
    Remove escape characters that colorize the text in a string.

    Parameters
    ----------
    text
        String to alter

    Returns
    -------
        String with all color escape characters removed
    """
    return re.compile(r'\x1b[^m]*m').sub("", text)


class StreamFormatter(logging.Formatter):
    """
    Custom `Formatter <logging.Formatter>` for the stream handler.
    """
    base_level = None
    """
    The base logging level for the class.  Used to determine whether or not to
    include the calling frame in the log message.
    """

    def format(self, record):

        # RGB colors for the logging levels
        level_colors = {
            'debug': [116, 173, 209],
            'info': [49, 54, 149],
            'test': [223, 255, 0],
            'warning': [253, 174, 97],
            'error': [215, 48, 39],
            'critical': [165, 0, 38],
        }
        frame_color = level_colors['debug']

        rec = copy.copy(record)
        levelname = rec.levelname.lower()
        if levelname not in level_colors:
            # Unknown level name so default to the standard formatter
            return logging.Formatter.format(self, record)

        # Add the level in colored text
        msg = color_text(f'[{levelname.upper()}]', level_colors[levelname], bold=True, nchar=10)
        msg += ' - '
        if self.base_level == logging.DEBUG:
            # If including debug messages, include file frame inspection in
            # *all* log messages.
            msg += color_text(f'{rec.filename}:{rec.funcName}:{rec.lineno}', frame_color) + ' - '
        # Add the message header
        rec.msg = msg + rec.msg

        # Return the base formatting
        return logging.Formatter.format(self, rec)


class DebugStreamFormatter(StreamFormatter):
    """
    Set the base logging level to DEBUG
    """
    base_level = logging.DEBUG


class FileFormatter(logging.Formatter):
    """
    Custom `Formatter <logging.Formatter>` for the file handler.
    """

    base_fmt = "%(levelname)8s | %(asctime)s | %(filename)s:%(funcName)s:%(lineno)s | %(message)s"

    def __init__(self, fmt=base_fmt):
        logging.Formatter.__init__(self, fmt, datefmt='%Y-%m-%d %H:%M:%S')


class PypeItLogger(logging.Logger):
    """
    Custom logging system for pypeit.

    This borrows heavily from implementations in astropy and sdsstools.
    """
    _excepthook_orig = None

    def init(self,
        level: int = logging.INFO,
        stream: io.TextIOBase | None = None,
        log_file: str | Path | None = None,
        log_file_level: int | None = None,
    ):
        """
        Initialise the logger.

        Parameters
        ----------
        level
            The logging level printed to the console
        stream
            Stream for logging messages, which defaults to sys.stderr.
        log_file
            Name for a log file.  If None, logging is only recorded to the
            console.  If the file provided already exists, it will be
            overwritten!
        log_file_level
            The logging level specific to the log file.  If None, adopt the
            console logging level.
        """
        # NOTE: I originally included these as options in the class.  I've
        # removed them for now (i.e., we'll always catch warnings and
        # exceptions), but I've left the if statements in place below in case we
        # want to make these things options in the future.
        capture_exceptions = True
        capture_warnings = True

        # NOTE: Because of how get_logger works, this makes warnings_logger an
        # instance of PypeItLogger.
        self.warnings_logger = logging.getLogger("py.warnings")

        # Set the base level of the logger to DEBUG    
        self.setLevel(logging.DEBUG)

        # Clear handlers before recreating.
        for handler in self.handlers.copy():
            if handler in self.warnings_logger.handlers:
                # Remove any added to the warnings logger
                self.warnings_logger.removeHandler(handler)
            self.removeHandler(handler)

        # Reset the exception hook (only if it was reset by this logger)
        if self._excepthook_orig is not None and sys.excepthook == self._excepthook:
            sys.excepthook = self._excepthook_orig
            self._excepthook_orig = None

        # Catch and parse exceptions
        if capture_exceptions:
            self._excepthook_orig = sys.excepthook
            sys.excepthook = self._excepthook

        # Set the stream handler, its formatting, its level, and then add it to
        # the set of handlers
        self.sh = logging.StreamHandler(stream=stream)
        formatter = DebugStreamFormatter() if level <= logging.DEBUG else StreamFormatter()
        self.sh.setFormatter(formatter)
        self.sh.setLevel(level)
        self.addHandler(self.sh)

        if capture_warnings:
            logging.captureWarnings(True)

            # Only enable the sh handler if none is attached to the warnings
            # logger yet. Prevents duplicated prints of the warnings.
            for handler in self.warnings_logger.handlers:
                if isinstance(handler, logging.StreamHandler):
                    return

            self.warnings_logger.addHandler(self.sh)

        # Get the file handler
        if log_file is None:
            self.fh = None
            self.log_filename = None
        else:
            if log_file_level is None:
                log_file_level = level
            self.log_file = Path(log_file).absolute()
            self.fh = logging.FileHandler(str(self.log_file), mode='w')
            self.fh.setFormatter(FileFormatter())
            self.fh.setLevel(log_file_level)
            self.addHandler(self.fh)

            if self.warnings_logger:
                self.warnings_logger.addHandler(self.fh)

    def _excepthook(self, etype, value, trace):
        """
        Override the default exception hook to log an error message.
        """
        tb = trace
        if tb is None:
            exc_info = None
        else:
            # If the traceback is available, jump to the calling frame, which
            # gets passed to makeRecord
            while tb.tb_next:
                tb = tb.tb_next
            exc_info = (etype, value, tb)

        # Add the error type to the message.
        if len(value.args) > 0:
            message = f"{etype.__name__}: {str(value)}"
        else:
            message = str(etype.__name__)

        # Log the error
        self.error(message, exc_info=exc_info)

        # Call the original exception hook
        self._excepthook_orig(etype, value, trace)

    @staticmethod        
    def convert_verbosity_to_logging_level(v):
        """
        Given a PypeIt "verbosity level," return the logging level.

        Parameters
        ----------
        v : int
            PypeIt verbosity level (0, 1, or 2)

        Returns
        -------
        int
            Corresponding logging level

        Raises
        ------
        ValueError
            Raised if the input verbosity level is not 0, 1, or 2.
        """
        match v:
            case 0:
                return logging.WARNING
            case 1:
                return logging.INFO
            case 2:
                return logging.DEBUG
            case _:
                raise ValueError(f'Verbosity level must be 0, 1, or 2, not {v}.')

    def makeRecord(
            self, name, level, pathname, lineno, msg, args, exc_info, func=None, extra=None,
            sinfo=None
    ):
        """
        Override the default makeRecord function to rework the message for exceptions.
        """

        # If the warning was issued by "warnings", try to recover the calling
        # frame details
        if name == 'py.warnings':
            frame = inspect.currentframe()
            save_frame = None
            while frame is not None:
                # Work backwards through the frame to find the first occurrence
                # of the call to the warnings.warn function.
                if (
                    Path(frame.f_code.co_filename).name == "warnings.py"
                    and frame.f_code.co_name == '_showwarnmsg'
                ):
                    save_frame = frame.f_back
                frame = frame.f_back
            if save_frame is not None:
                pathname = save_frame.f_code.co_filename
                lineno = save_frame.f_lineno
                func = save_frame.f_code.co_name

        # Do the same if (1) this is an error message, (2) the execution
        # information is provided, and (3) the error originates from the
        # exception hook.
        elif (level == logging.ERROR
            and exc_info is not None
            and Path(pathname).name == 'logger.py'
            and func is not None
            and func == '_excepthook'
        ):
            calling_frame = traceback.extract_tb(exc_info[2])[-1]
            pathname = calling_frame.filename
            lineno = calling_frame.lineno
            func = calling_frame.name
            # This keeps the traceback from being printed twice!
            exc_info = None

        # Call the base-class method
        return logging.Logger.makeRecord(
            self, name, level, pathname, lineno, msg, args, exc_info, func=func, extra=extra,
            sinfo=sinfo
        )
    
    def close_file(self):
        """
        Explicitly close the log file.
        """
        if self.fh is None:
            return
        self.fh.close()
        self.removeHandler(self.fh)
        if self.fh in self.warnings_logger.handlers:
            self.warnings_logger.removeHandler(self.fh)


# NOTE: If we allow warning and exception capture to be optional, remember to
# add them as parameters here as well.
def get_logger(
    level: int = logging.INFO,
    stream: io.TextIOBase | None = None,
    log_file: str | Path | None = None,
    log_file_level: int | None = None,
) -> PypeItLogger:
    """
    Instantiate a new logger.

    Parameters
    ----------
    level
        The logging level printed to the console
    stream
        Stream for logging messages, which defaults to sys.stderr.
    log_file
        Name for a log file.  If None, logging is only recorded to the
        console.  If the file provided already exists, it will be
        ovewritten!
    log_file_level
        The logging level specific to the log file.  If None, adopt the
        console logging level.

    Returns
    -------
        Logging object for PypeIt.
    """

    orig_logger = logging.getLoggerClass()
    logging.setLoggerClass(PypeItLogger)

    try:
        log = logging.getLogger("pypeit")
        log.init(
            level=level,
            stream=stream,
            log_file=log_file,
            log_file_level=log_file_level
        )
        # TODO: We might want to prohibit propagation of this logger to the root
        # one, but I'm not really sure if that's necessary or how it works.
        # Leaving this commented out for now.
        # log.propagate = False
    finally:
        logging.setLoggerClass(orig_logger)

    return log
