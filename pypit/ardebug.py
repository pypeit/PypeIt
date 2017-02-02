from __future__ import (print_function, absolute_import, division, unicode_literals)

from pypit import ginga as pyp_g

def init():
    """
    Returns
    -------
    debug : dict
        default debug dict
    """
    debug = dict(develop=False,
                 arc=False,
                 obj_profile=False,
                 sky_sub=False,
                 trace=False,
                 wave=False,
                 tilts=False,
                 flexure=False,
                 no_qa=False,
                 trace_obj=False,
                 )
    return debug


def show_image(args, **kwargs):
    """ Wrapper to pypit_ginga
    Parameters
    ----------
    args
    kwargs

    Returns
    -------

    """
    return pyp_g.show_image(args, **kwargs)

from pdb import *
