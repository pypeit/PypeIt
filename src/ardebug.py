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
                 new_instrument=False,
                 arc=False,
                 trace=False,
                 wave=False,
                 trace_obj=False,
                 )
    return debug
