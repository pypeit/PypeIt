""" Some PYPIT utilities for the package
"""
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

import re
from os.path import realpath, join

this_file = realpath(__file__)
this_path = this_file[:this_file.rfind('/')]


def get_version():
    """Get the value of ``__version__`` without having to import the module.
    Copied from DESI desiutils

    Parameters
    ----------

    Returns
    -------
    :class:`str`
        The value of ``__version__``.
    """
    ver = 'unknown'
    version_file = join(this_path, '_version.py')
    with open(version_file, "r") as f:
        for line in f.readlines():
            mo = re.match("__version__ = '(.*)'", line)
            lu = re.match("__lastupdate__ = '(.*)'", line)
            if mo:
                ver = mo.group(1)
            if lu:
                upd = lu.group(1)
    return ver, upd


def get_dummy_logger():
    """ Useful for testing

    Returns
    -------

    """
    import ardebug
    import armsgs as pyparm
    debug = ardebug.init()

    version, last_updated = get_version()
    pyparm.pypit_logger = pyparm.Messages(None, debug, last_updated, version, 1)
    return pyparm.pypit_logger
