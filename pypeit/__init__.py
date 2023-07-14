"""
pypeit package initialization.

The current main purpose of this is to provide package-level globals
that can be imported by submodules.
"""

# Imports for signal and log handling
import os
import sys
import signal
import warnings

from .version import version

def short_warning(message, category, filename, lineno, file=None, line=None):
    """
    Return the format for a short warning message.
    """
    return ' %s: %s (%s:%s)\n' % (category.__name__, message, os.path.split(filename)[1], lineno)

warnings.formatwarning = short_warning


# Set version
__version__ = version

# Report current coverage
__coverage__ = 0.55

# Import and instantiate the logger
from pypeit import pypmsgs
msgs = pypmsgs.Messages()

# Import the close_qa method so that it can be called when a hard stop
# is requested by the user
from pypeit.core.qa import close_qa


# Send all signals to messages to be dealt with (i.e. someone hits ctrl+c)
def signal_handler(signalnum, handler):
    """
    Handle signals sent by the keyboard during code execution
    """
    if signalnum == 2:
        msgs.info('Ctrl+C was pressed. Ending processes...')
        close_qa(msgs.pypeit_file, msgs.qa_path)
        msgs.close()
        sys.exit()


signal.signal(signal.SIGINT, signal_handler)

# Ignore all warnings given by python
# TODO: I'd rather we not do this.  Is there a way we can redirect
# warnings to pypeit.msgs ?
#warnings.resetwarnings()
#warnings.simplefilter('ignore')


