""" Module for ginga routines.  Mainly for debugging
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import numpy as np

from pypit import armsgs

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


def connect_to_ginga(host='localhost', port=9000,):
    """ Connect to an active RC Ginga
    Parameters
    ----------
    host : str, optional
    port : int, optional

    Returns
    -------
    viewer : RemoteClient
      connectoin to Ginga

    """
    from ginga.util import grc as ggrc
    # Start
    viewer = ggrc.RemoteClient(host, port)
    # Return
    return viewer


def show_image(img):
    """ Displays input image in Ginga viewer
    Supersedes method in xastropy

    Parameters
    ----------
    img : ndarray (2D)

    Returns
    -------

    """
    viewer = connect_to_ginga()
    ch = viewer.channel('Image')
    name='image'
    ch.load_np(name, img, 'fits', {})


def chk_arc_tilts(msarc, trcdict, sedges=None, yoff=0., xoff=0.):
    """  Display arc image and overlay the arcline tilt measurements
    Parameters
    ----------
    msarc : ndarray
    trcdict : dict
      Contains trace info
    sedges : tuple
      Arrays of the slit
    xoff : float, optional
      In case Ginga has an index offset.  It appears not to
    yoff : float, optional


    Returns
    -------

    """
    # Connect
    cname = 'ArcTilts'
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
            points = list(zip(x[gdy].tolist(),y[gdy].tolist()))
            if trcdict['aduse'][idx]:
                clr = 'green'  # Good line
            else:
                clr = 'red'  # Bad line
            canvas.add('path', points, color=clr)
    msgs.info("Check the Ginga viewer")
    # Show slit edges
    if sedges is not None:
        y = (np.arange(msarc.shape[0]) + yoff).tolist()
        # Left
        for edge in [0,1]:
            points = zip(sedges[edge].tolist(),y)
            canvas.add('path', points, color='cyan')
    # ALTERNATE for comparing methods
    if 'save_yt' in trcdict.keys():
        x = trcdict['xtfit'][2]
        y = trcdict['save_yt']
        gdy = y > 0.
        points = list(zip(x[gdy].tolist(),y[gdy].tolist()))
        canvas.add('path', points, color='blue')
    #debugger.set_trace()

