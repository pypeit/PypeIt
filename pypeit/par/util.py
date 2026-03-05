# -*- coding: utf-8 -*-
"""
Utility functions for PypeIt parameter sets

.. include:: ../include/links.rst
"""
import ast
import itertools
from typing import Type

from configobj import ConfigObj
from IPython import embed

from pypeit import PypeItError


def _eval_ignore():
    """Provides a list of strings that should not be evaluated."""
    return [ 'open', 'file', 'dict', 'list', 'tuple' ]


def ast_literal_eval(inp):
    """
    A wrapper for :func:`ast.literal_eval` that returns the input if it raises a
    ValueError.
    """
    try:
        return ast.literal_eval(inp)
    except ValueError:
        return inp


def _eval_iter(inp:list[str], left:str, right:str, otype:Type) -> list:
    """
    Convenience function used to abstract the core functionality of
    :func:`eval_tuple` and :func:`eval_list`.

    Parameters
    ----------
    inp
        Input list of strings
    left
        Left bracket delimeter
    right
        Right bracket delimeter
    otype
        Return type for the components of the list

    Returns
    -------
        A list of objects with type ``otype`` as converted from the provided
        strings.
    """
    grps = [
        i for i in list(itertools.chain(*[s.split(right) for s in ','.join(inp).split(left)]))
        if len(i) > 1
    ]
    return list(otype(map(ast_literal_eval, g.split(','))) for g in grps)


def eval_tuple(inp:list[str]) -> list[tuple]:
    """
    Evaluate the input to one or more tuples.

    This allows conversion of one or more tuples provided to a configuration
    parameters.

    .. warning::

        - Currently can only handle tuples with integers, floats, or strings!

    Parameters
    ----------
    inp 
        A list of strings that are converted into a list of tuples.  The
        parentheses must be within the list of elements.

    Returns
    -------
        A list of tuples with the converted elements.
    """
    return _eval_iter(inp, '(', ')', tuple)


def eval_list(inp:list[str]) -> list[list]:
    """
    Evaluate the input to one or more lists.

    This allows conversion of one or more lists provided to a configuration
    parameters.

    .. warning::

        - Currently can only handle lists with integers, floats, or strings!

    Parameters
    ----------
    inp 
        A list of strings that are converted into a list of lists.  The square
        brackets must be within the list of elements.

    Returns
    -------
        A list of lists with the converted elements.
    """
    return _eval_iter(inp, '[', ']', list)


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
                d[k] = eval_tuple(d[k])
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
                        replacement += [ ast_literal_eval(v) ]
                    except:
                        replacement += [ v ]
            d[k] = replacement
            continue

        try:
            d[k] = ast_literal_eval(d[k]) if d[k] not in ignore else d[k]
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
