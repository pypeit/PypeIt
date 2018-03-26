"""
pypit package initialization.

The current main purpose of this is to provide package-level globals
that can be imported by submodules.
"""

# Set version
__version__ = '0.8.0.dev0'
# TODO: Do we need this?  A release should be tied to a date, and the
# date information could be put in the CHANGES.rst
__last_updated__ = '23Mar2018'

# Imports for singal and log handling
import sys
import signal
import warnings

# Import and instantiate the logger
from pypit import armsgs
msgs = armsgs.Messages()

# Import the close_qa method so that it can be called when a hard stop
# is requested by the user
# TODO: Is this necessary?  As far as I can tell, close_qa() doesn't
# actually close any already opened qa files
from pypit.arqa import close_qa

# Send all signals to messages to be dealt with (i.e. someone hits ctrl+c)
def signal_handler(signalnum, handler):
    """
    Handle signals sent by the keyboard during code execution
    """
    if signalnum == 2:
        msgs.info('Ctrl+C was pressed. Ending processes...')
        close_qa(msgs.pypit_file)
        msgs.close()
        sys.exit()

signal.signal(signal.SIGINT, signal_handler)

# Ignore all warnings given by python
warnings.resetwarnings()
warnings.simplefilter('ignore')

