""" Module for ginga routines.  Mainly for debugging
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import numpy as np

from pypit import armsgs
from pypit import pyputils

# CANNOT LOAD DEBUGGER AS THIS MODULE IS CALLED BY ARDEBUG
#from pypit import ardebug as debugger
import pdb as debugger

try:
    basestring
except NameError:
    basestring = str

# Logging
msgs = armsgs.get_logger()   # THESE MAY NOT WORK..
'''
if msgs is None:  # For usage outside of PYPIT
    msgs = pyputils.get_dummy_logger()
    armsgs.pypit_logger = msgs
'''


def connect_to_ginga(host='localhost', port=9000):
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
    # Test
    ginga = viewer.shell()
    try:
        tmp = ginga.get_current_workspace()
    except:
        msgs.warn("Problem connecting to Ginga.  Launch an RC Ginga viewer then continue.")
        debugger.set_trace()
    # Return
    return viewer


def show_image(inp, chname='Image', **kwargs):
    """ Displays input image in Ginga viewer
    Supersedes method in xastropy

    Parameters
    ----------
    inp : str or ndarray (2D)
      If str, assumes the image is written to disk

    Returns
    -------

    """
    from astropy.io import fits
    from pypit import arlris
    if isinstance(inp, basestring):
        if '.fits' in inp:
            if 'raw_lris' in kwargs.keys():
                img, head, _ = arlris.read_lris(inp)
            else:
                hdu = fits.open(inp)
                try:
                    exten = kwargs['exten']
                except KeyError:
                    exten = 0
                img = hdu[exten].data
    else:
        img = inp

    viewer = connect_to_ginga()
    ch = viewer.channel(chname)
    name='image'
    # Header
    header = {}
    header['NAXIS1'] = img.shape[1]
    header['NAXIS2'] = img.shape[0]
    header['WCS-XIMG'] = '/home/xavier/REDUX/Keck/LRIS/2017mar20/lris_red_setup_C/MF_lris_red/MasterWave_C_02_aa.fits'
    # Giddy up
    ch.load_np(name, img, 'fits', header)
    return viewer, ch


def show_slits(viewer, ch, lordloc, rordloc, slit_ids):
    """ Overplot slits on image in Ginga
    Parameters
    ----------
    viewer
    ch
    lordloc
    rordloc
    slit_ids : list of int

    Returns
    -------

    """
    # Canvas
    canvas = viewer.canvas(ch._chname)
    canvas.clear()
    # y-axis
    y = (np.arange(lordloc.shape[0])).tolist()
    #ohf = lordloc.shape[0] // 2
    tthrd = int(2*lordloc.shape[0]/3.)
    # Loop on slits
    for slit in range(lordloc.shape[1]):
        # Left
        points = list(zip(lordloc[:,slit].tolist(),y))
        canvas.add('path', points, color='cyan')
        # Right
        points = list(zip(rordloc[:,slit].tolist(),y))
        canvas.add('path', points, color='cyan')
        # Text -- Should use the 'real' name
        canvas.add('text', float(lordloc[tthrd,slit]), float(y[tthrd]),
                   'S{:d}'.format(slit_ids[slit]), color='cyan')

def show_trace(viewer, ch, trace, trc_name, color='blue', clear=False):
    # Canvas
    canvas = viewer.canvas(ch._chname)
    if clear:
        canvas.clear()
    # Show
    y = (np.arange(trace.size)).tolist()
    points = list(zip(trace.tolist(),y))
    canvas.add('path', points, color=color)
    # Text
    ohf = trace.size // 2
    canvas.add('text', float(trace[ohf]), float(y[ohf]), trc_name,
               rot_deg=90., color=color, fontsize=17.)

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
    canvas = viewer.canvas(ch._chname)
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
    #msgs.info("Check the Ginga viewer")
    # Show slit edges
    if sedges is not None:
        y = (np.arange(msarc.shape[0]) + yoff).tolist()
        # Left
        for edge in [0,1]:
            points = list(zip(sedges[edge].tolist(),y))
            canvas.add('path', points, color='cyan')
    # ALTERNATE for comparing methods
    if 'save_yt' in trcdict.keys():
        x = trcdict['xtfit'][2]
        y = trcdict['save_yt']
        gdy = y > 0.
        points = list(zip(x[gdy].tolist(),y[gdy].tolist()))
        canvas.add('path', points, color='blue')
    #debugger.set_trace()

