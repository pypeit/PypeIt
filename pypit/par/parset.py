# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Define a utility base class used to hold parameters.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Class usage examples*:
    to be added

.. todo::
    - Add range and length parameters allowing one to define the range
      allowed for the parameter values and number of elements required
      (if the parameter is an array)
    - Allow for a from_par_file classmethod to initialize the parameter
      set based on a yanny parameter file.
    - Save the defaults and allow for a revert_to_default function.
    - Write an __add__ function that will all you to add multiple
      parameter sets.

*Revision history*:
    | **16 Jun 2015**: Original implementation by K. Westfall (KBW)
    | **18 Mar 2016**: (KBW) Change dtype checking
    | **23 Mar 2016**: (KBW) Changed initialization type checking of
        lists to use `isinstance`_.
    | **02 Apr 2016**: (KBW) Allow input parameters to be callable
        functions.
    | **05 Apr 2018**: (KBW) Added to pypit repo

.. _isinstance: https://docs.python.org/2/library/functions.html#isinstance

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy

class ParSet(object):
    """
    Generic base class to handle and manipulate a list of operational
    parameters.  A glorified dictionary that constrains and types its
    components.

    Args:
        pars (list) : A list of keywords for a list of parameter values
        values (list) : (**Optional**) Initialize the parameters to
            these values.  If not provided, all parameters are
            initialized to `None` or the provided default.
        defaults (list) : (**Optional**) For any parameters not provided
            in the *values* list, use these default values.  If not
            provided, no defaults are assumed.
        options (list) : (**Optional**) Force the parameters to be one
            of a list of options.  Each element in the list can be a
            list itself.  If not provided, all parameters are allowed to
            take on any value within the allowed data type.
        dtypes (list) : (**Optional**) Force the parameter to be one of
            a list of data types.  Each element in the list can be a
            list itself.  If not provided, all parameters are allowed to
            have any data type.
        can_call (list) : (**Optional**) Flag that the parameters are
            callable operations.  Default is False.

    Raises:
        TypeError: Raised if the input parameters are not lists or if
            the input keys are not strings.
        ValueError: Raised if any of the optional arguments do not have
            the same length as the input list of parameter keys.

    Attributes:
        npar (int) : Number of parameters
        data (dict) : Dictionary with the parameter values
        options (dict) : Dictionary with the allowed options for the
            parameter values
        dtype (dict) : Dictionary with the allowed data types for the
            parameters
        can_call (dict): Dictionary with the callable flags

    """
    def __init__(self, pars, values=None, defaults=None, options=None, dtypes=None, can_call=None):
        # Check that the list of input parameters is a list of strings
        if not isinstance(pars, list):
            raise TypeError('Input parameter keys must be provided as a list.')
        for key in pars:
            if not isinstance(key, str):
                raise TypeError('Input parameter keys must be strings.')
        
        # Get the length of the parameter list and make sure the list
        # has unique values
        self.npar = len(pars)
        if len(numpy.unique(numpy.array(pars))) != self.npar:
            raise ValueError('All input parameter keys must be unique.')

        # Check that the other lists, if provided, have the correct type
        # and length
        if values is not None and (not isinstance(values, list) or len(values) != self.npar):
            raise ValueError('Values must be a list with the same length as the keys list.')
        if defaults is not None and (not isinstance(defaults, list) or len(defaults) != self.npar):
            raise ValueError('Defaults must be a list with the same length as the keys list.')
        if options is not None and (not isinstance(options, list) or len(options) != self.npar):
            raise ValueError('Options must be a list with the same length as the keys list.')
        if dtypes is not None and (not isinstance(dtypes, list) or len(dtypes) != self.npar):
            raise ValueError('Data types list must have the same length as the keys list.')
        if can_call is not None and (not isinstance(can_call, list) or len(can_call) != self.npar):
            raise ValueError('List of callable flags must have the same length as keys list.')

        # Set up dummy lists for no input
        if values is None:
            values = [None]*self.npar
        if defaults is None:
            defaults = [None]*self.npar
        if options is None:
            options = [None]*self.npar
        if dtypes is None:
            dtypes = [None]*self.npar
        if can_call is None:
            can_call = [False]*self.npar

        # Set the valid options
        self.options = dict([ (p, [o]) if o is not None and not isinstance(o, list) else (p, o) \
                                       for p, o in zip(pars, options) ])
        # Set the valid types
        self.dtype = dict([ (p, [t]) if t is not None and not isinstance(t, list) else (p, t) \
                                     for p, t in zip(pars, dtypes) ])

        # Set the calling flags
        self.can_call = dict([ (p, t) for p, t in zip(pars, can_call) ])

        # Set the data dictionary using the internal functions
        self.data = {}
        for p, d, v in zip(pars, defaults, values):
            if v is None:
                self.__setitem__(p, d)
                continue
            self.__setitem__(p, v)


    def __getitem__(self, key):
        """
        Return the value of the designated key.

        Args:
            key (str) : Key for new parameter
        """
        return self.data[key]


    def __setitem__(self, key, value):
        """
        Set the value for a key.

        Args:
            key (str) : Key for new parameter
            value (*dtype*) : Parameter value, must have a type in the
                list provided (*dtype*), if the list is provided

        Raises:
            ValueError: Raised if the parameter value is not among the
                allowed options (:attr:`options`).
            TypeError: Raised if the parameter value does not have an
                allowed data type (:attr:`dtype`) or if the provided
                value is not a callable object, but is expected to be by
                :attr:`can_call`.
        """
        if value is None:
            self.data[key] = value
            return

        if self.options[key] is not None and value not in self.options[key]:
            raise ValueError('Input value invalid: {0}.\nOptions are: {1}'.format(value,
                                                                                self.options[key]))
        if self.dtype[key] is not None \
                and not any([ isinstance(value, d) for d in self.dtype[key]]):
            raise TypeError('Input value incorrect type: {0}.\nValid types are: {1}'.format(value,
                                                                                self.dtype[key]))

        if self.can_call[key] and not callable(value):
            raise TypeError('{0} is not a callable object.'.format(value))

        self.data[key] = value


    def __len__(self):
        """Return the number of parameters."""
        return self.npar
        

    def __iter__(self):
        """Return an iterable to the parameter values."""
        return iter(self.data.values())


    def __repr__(self):
        """Return a crude string represenation of the parameters."""
        out = ''
        for k in self.data.keys():
            out += '{0:>10}: {1}\n'.format(k, self.data[k])
        return out


    def keys(self):
        return list(self.data.keys())

    
    def add(self, key, value, options=None, dtype=None, can_call=None):
        """
        Add a new parameter.

        Args:
            key (str) : Key for new parameter
            value (*dtype*) : Parameter value, must have a type in the
                list provided (*dtype*), if the list is provided
            options (list) : (**Optional**) List of discrete values that
                the parameter is allowed to have
            dtype (list) : (**Optional**) List of allowed data types
                that the parameter can have
            can_call (bool) : (**Optional**) Flag that the parameters
                are callable operations.  Default is False.
        """
        if key in self.data.keys():
            raise ValueError('Keyword {0} already exists and cannot be added!')
        self.npar += 1
        self.options[key] = [options] if options is not None and not isinstance(options, list) \
                                      else options
        self.dtype[key] = [dtype] if dtype is not None and not isinstance(dtype, list) else dtype
        self.can_call[key] = False if can_call is None else can_call
        try:
            self.__setitem__(key, value)
        except:
            # Delete the added components
            del self.options[key]
            del self.dtype[key]
            del self.can_call[key]
            # Re-raise the exception
            raise


class ParDatabase(object):
    """

    Class used as a list of ParSets in a glorified structured numpy
    array.

    Very similar to yanny when converted to a numpy array.

    Can be initialized using a list of ParSet objects, or an SDSS
    parameter file.

    .. todo::

        - Check that the data types are the same for all ParSet objects
          in the list
        - Better handle the NaN values when converting None to a float
          type
        - Add from_par_file classmethod?
    
    """
    def __init__(self, inp):
        """
        nsets - number of parameter sets
        npar - number of parameters in each parameter set
        data - parameter set data
        options - allowed options for values
        dtype - allowed datatypes
        can_call - parameter is a callable function
        """
        _inp = [inp] if isinstance(inp, ParSet) else inp
        if not isinstance(_inp, list):
            raise TypeError('Input must be a list.')
        for i in _inp:
            if not isinstance(i, ParSet):
                raise TypeError('Input must be a list of ParSet objects.')
        self.npar = _inp[0].npar
        self.nsets = len(_inp)
        keys = _inp[0].keys()
        for i in range(1,self.nsets):
            if _inp[i].npar != self.npar:
                raise ValueError('Not all ParSet objects have the same number of parameters.')
            if _inp[i].keys() != keys:
                raise ValueError('Not all ParSet objects have the same keys.')
            # Other checks?

        record_dtypes = self._set_dtypes(_inp, 0)

        data = []
        for i in range(self.nsets):
            data += [ tuple([_inp[i][k] for k in keys]) ]

        # WARNING: None values are converted to nan if data type is
        # float
        self.data = numpy.array(data, dtype=record_dtypes ).view(numpy.recarray)
        self.options = inp[0].options.copy()
        self.dtype = inp[0].dtype.copy()
        self.can_call = inp[0].can_call.copy()
   

    def __getitem__(self, key):
        """
        Return the value of the designated key.

        Args:
            key (str) : Key for new parameter
        """
        return self.data[key]


    @staticmethod
    def _set_dtypes(inp, i):
        keys = inp[i].keys()
        dtypes = []
        for k in keys:
            if inp[i].dtype[k] is None:
                dtypes += [(k,object)]
                continue
            # inp.dtype is always a list
            if any([t in inp[i].dtype[k] for t in [int , float]]) \
                and any([t in inp[i].dtype[k] for t in [list, numpy.ndarray]]):
                warnings.warn('Parameter set has elements that can be either individual ' \
                              'ints/floats or lists/arrays.  Database column {0} will have type ' \
                              '\'object\'.'.format(k))
                dtypes += [(k,object)]
            elif len(list({int, float} - set(inp[i].dtype[k]))) == 0:
                dtypes += [(k,float)]
            elif len(list({list, numpy.ndarray} - set(inp[i].dtype[k]))) == 0 \
                    or inp[i].dtype[k] == numpy.ndarray:
                _inp = numpy.asarray(inp[i][k])
                dtypes += [(k,_inp.dtype,_inp.shape)]
            elif isinstance(inp[i][k], str):
                if any([ _inp[k] is None for _inp in inp]):
                    dtypes += [(k, object)]
                else:
                    dtypes += [(k,'<U{0:d}'.format(max([ len(_inp[k]) for _inp in inp])))]
            else:
                dtypes += [(k,type(inp[i][k]))]
        return dtypes


    def append(self, pdb):
        if not isinstance(pdb, ParDatabase):
            raise TypeError('Can only append ParDatabase object.')

        try:
            self.data = numpy.append(self.data, pdb.data)
        except TypeError as e:
            raise TypeError('Could not append data:: {0}'.format(e))
            



