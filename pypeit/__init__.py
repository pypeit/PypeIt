"""
pypeit package initialization.

The current main purpose of this is to provide package-level globals
that can be imported by submodules.
"""

from .pkginit.version import version

# Set version
__version__ = version

# Report current coverage
# TODO: How old is this?  Can we update it automatically?
__coverage__ = 0.55

# Start the log
import logging
from .pkginit.logger import get_logger
log = get_logger(level=logging.DEBUG)

# Import and instantiate the data path parser
# NOTE: This *MUST* come after log and __version__ are defined above
#from pypeit import pypeitdata
#from .pkginit.pypeitdata import PypeItDataPaths
from .pkginit.cache import PypeItDataPaths
dataPaths = PypeItDataPaths()

# Expose all the exceptions
from .pkginit.exceptions import *

## Imports for signal and log handling
#import sys
#import signal
## Send all signals to messages to be dealt with (i.e. someone hits ctrl+c)
#def signal_handler(signalnum, handler):
#    """
#    Handle signals sent by the keyboard during code execution
#    """
#    if signalnum == 2:
#        log.info('Ctrl+C was pressed. Ending processes...')
#        sys.exit()
#
#signal.signal(signal.SIGINT, signal_handler)

