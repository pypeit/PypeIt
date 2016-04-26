import sys
from os.path import dirname, basename
from textwrap import wrap as wraptext
from inspect import currentframe, getouterframes
from glob import glob

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
                 trace_obj=False,
                 wave=False,
                 )
    return debug
