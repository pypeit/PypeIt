"""
pypeit package initialization.

The current main purpose of this is to provide package-level globals
that can be imported by submodules.
"""

from .pkg.version import version

# Set version
__version__ = version

# Report current coverage
# TODO: How old is this?  Can we update it automatically?
__coverage__ = 0.55

# Start the log
import logging
from .pkg.logger import get_logger
log = get_logger(level=logging.DEBUG)

# Import and instantiate the data path parser
# NOTE: This *MUST* come after log and __version__ are defined above
from .pkg.pypeitdata import PypeItDataPaths
dataPaths = PypeItDataPaths()

# Import all the exceptions so that they can be directly imported (e.g., `from
# pypeit import PypeItError`) in all package imports.
from .pkg.exceptions import *

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

