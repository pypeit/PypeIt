from astropy.io import fits
from astropy.units import Quantity

import armsgs

# Logging
msgs = armsgs.get_logger()

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

