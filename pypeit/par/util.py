# -*- coding: utf-8 -*-
"""
Utility functions for PypeIt parameter sets

.. include:: ../include/links.rst
"""

from configobj import ConfigObj
from IPython import embed

from pypeit import PypeItError
from pypeit import utils


def _eval_ignore():
    """Provides a list of strings that should not be evaluated."""
    return [ 'open', 'file', 'dict', 'list', 'tuple' ]


def recursive_dict_evaluate(d):
    """
    Recursively run :func:`eval` on each element of the provided
    dictionary.

    A raw read of a configuration file with `configobj`_ results in a
    dictionary that contains strings or lists of strings.  However, when
    assigning the values for the various ParSets, the `from_dict`
    methods expect the dictionary values to have the appropriate type.
    E.g., the `configobj`_ will have something like d['foo'] = '1', when
    the `from_dict` method expects the value to be an integer (d['foo']
    = 1).

    This function tries to evaluate *all* dictionary values, except for
    those listed above in the :func:`_eval_ignore` function.  Any value
    in this list or where::

        ast_literal_eval(d[k]) for k in d.keys()

    raises an exception is returned as the original string.

    This is currently only used in
    :func:`~pypeit.par.pypeitpar.PypitPar.from_cfg_file`; see further comments
    there.

    Args:
        d (dict):
            Dictionary of values to evaluate

    Returns:
        dict: Identical to input dictionary, but with all string values
        replaced with the result of `eval(d[k])` for all `k` in
        `d.keys()`.
    """
    ignore = _eval_ignore()
    for k in d.keys():
        if isinstance(d[k], dict):
            # Recursive call to deal with nested dictionaries
            d[k] = recursive_dict_evaluate(d[k])
            continue

        if isinstance(d[k], list) and any(['(' in e for e in d[k]]):
            # NOTE: This enables syntax for constructing one or more tuples.
            try:
                d[k] = utils.eval_tuple(d[k])
            except (PypeItError, SyntaxError) as e:
                # The tuple evaluation failed.  Assume that this can be handled
                # later in the code and leave the dictionary element unaltered.
                # 
                # SyntaxError is raised for entries that include a tuple for the
                # mosaic and a series of locations in the mosaiced image, like
                # add_slits, rm_slits, and manual.
                pass
            continue

        if isinstance(d[k], list):
            replacement = []
            for v in d[k]:
                if v in ignore:
                    replacement += [ v ]
                else:
                    try:
                        replacement += [ utils.ast_literal_eval(v) ]
                    except:
                        replacement += [ v ]
            d[k] = replacement
            continue

        try:
            d[k] = utils.ast_literal_eval(d[k]) if d[k] not in ignore else d[k]
        except:
            pass

    return d


def parset_to_dict(par):
    """
    Convert the provided parset into a dictionary.

    Args:
        par (ParSet):

    Returns:
        dict: Converted ParSet

    """
    try:
        d = dict(ConfigObj(par.to_config(section_name='tmp'))['tmp'])
    except:
        d = dict(ConfigObj(par.to_config()))
    return recursive_dict_evaluate(d)
