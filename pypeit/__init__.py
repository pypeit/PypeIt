"""
pypeit package initialization.

The current main purpose of this is to provide package-level globals
that can be imported by submodules.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# Imports for signal and log handling
import sys
import signal
import warnings

# Set version
__version__ = '0.9.1dev'

# Import and instantiate the logger
from pypeit import pypmsgs
msgs = pypmsgs.Messages()

from pypeit import check_requirements  # THIS IMPORT DOES THE CHECKING.  KEEP IT

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
        close_qa(msgs.pypit_file)
        msgs.close()
        sys.exit()

signal.signal(signal.SIGINT, signal_handler)

# Ignore all warnings given by python
warnings.resetwarnings()
warnings.simplefilter('ignore')

