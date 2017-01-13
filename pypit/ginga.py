from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import re
import sys
import shutil
import string
import numpy as np
import yaml

from pypit import armsgs
from pypit import arparse as settings
from pypit import arutils
from pypit.arflux import find_standard_file
from astropy.io.votable.tree import VOTableFile, Resource, Table, Field
from astropy.table import Table as tTable, Column
from astropy import units as u

from linetools import utils as ltu

try:
    from xastropy.xutils import xdebug as debugger
except ImportError:
    import pdb as debugger
try:
    basestring
except NameError:
    basestring = str

# Logging
msgs = armsgs.get_logger()


def connect_to_ginga():
    from ginga.util import grc as ggrc
    host='localhost'
    port=9000
    # Start
    viewer = ggrc.RemoteClient(host, port)
    # Return
    return viewer


def show_image(img):
    viewer = connect_to_ginga()
    ch = viewer.channel('Image')
    name='image'
    ch.load_np(name, img, 'fits', {})


def chk_arc_tilts(msarc, trcdict, sedges=None, yoff=0., xoff=0.):
    # Connect
    cname = 'ATilts'
    viewer = connect_to_ginga()
    ch = viewer.channel(cname)
    canvas = viewer.canvas(cname)
    # Show image, clear canvas [in case this is a repeat]
    name='image'
    ch.load_np(name, msarc, 'fits', {})
    canvas.clear()
    # Show a trace
    ntrc = len(trcdict['arcdet'])
    for idx in range(ntrc):
        if trcdict['xtfit'][idx] is None:
            continue
        x = trcdict['xtfit'][idx] + xoff
        y = trcdict['ytfit'][idx] + yoff  # FOR IMAGING (ALREADY OFFSET IN GINGA)
        gdy = y > 0.
        if np.sum(gdy) > 0:
            points = zip(x[gdy].tolist(),y[gdy].tolist())
            if trcdict['aduse'][idx]:
                clr = 'green'
            else:
                clr = 'red'
            canvas.add('path', points, color=clr)
    msgs.info("Check the Ginga viewer")
    # Show slit edges
    if sedges is not None:
        y = (np.arange(msarc.shape[0]) + yoff).tolist()
        # Left
        for edge in [0,1]:
            points = zip(sedges[edge].tolist(),y)
            canvas.add('path', points, color='cyan')
    # 217.10 -> 217.49
    dy = trcdict['ytfit'][2] - np.roll(trcdict['ytfit'][2],1)
    bad = dy > 0.1
    print('dy', dy[bad])
    # ALTERNATE
    if 'save_yt' in trcdict.keys():
        x = trcdict['xtfit'][2]
        y = trcdict['save_yt']
        gdy = y > 0.
        points = zip(x[gdy].tolist(),y[gdy].tolist())
        canvas.add('path', points, color='blue')
    debugger.set_trace()

