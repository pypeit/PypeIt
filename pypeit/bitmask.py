# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Base class for handling bit masks.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Class usage examples*:
    TODO

*Revision history*:
    | **01 Jun 2015**: Original implementation by K. Westfall (KBW)
    | **07 Oct 2015**: (KBW) Added a usage case
    | **29 Jan 2016**: (KBW) Changed attributes of :class:`BitMask` and
        added functionality to print a description of the bits.  Convert
        :class:`mangadap.proc.templatelibrary.TemplateLibraryBitMask` to
        new format where the bits are read from a configuration file.
    | **17 Feb 2016**: (KBW) Minor edit to documentation
    | **16 Mar 2016**: (KBW) Moved TemplateLibraryBitMask to
        :class:`mangadap.proc.templatelibrary.TemplateLibraryBitMask`;
        moved HDUList_mask_wavelengths to
        :func:`mangadap.proc.util.HDUList_mask_wavelengths`.
    | **27 Mar 2016**: (KBW) Added :func:`BitMask.from_ini_file` and
        :func:`BitMask.from_par_file` class methods.  Added
        :func:`BitMask._fill_sequence` static method.  This allows for
        :class:`BitMask` objects to be declared directly from the files,
        and allows the bit values to take on any number.
    | **05 Apr 2016**: (KBW) Added parameters to initialization of
        :class:`BitMask` objects to clean up some of the derived class
        initialization.
    | **11 May 2016**: (KBW) Switch to using `pydl.pydlutils.yanny`_
        instead of internal yanny reader
    | **29 Jul 2016**: (KBW) Change asarray to atleast_1d
    | **06 Sep 2018**: (KBW) Added to PypeIt repo, removed
        functionality of instantiating a bitmask from a file, code
        update, and slight doc changes.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy
import os
import textwrap

class BitMask:
    r"""
    Generic class to handle and manipulate bitmasks.  The input list of
    bit names (keys) must be unique, except that values of 'NULL' are
    ignored.  The index in the input keys determines the bit value;
    'NULL' keys are included in the count.  For example::

        >>> from mangadap.util.bitmask import BitMask
        >>> keys = [ 'key1', 'key2', 'NULL', 'NULL', 'key3' ]
        >>> bm = BitMask(keys)
        >>> bm.info()
                 Bit: key1 = 0

                 Bit: key2 = 1

                 Bit: key3 = 4

    Args:
        keys (:obj:`str`, :obj:`list`):
            List of keys (or single key) to use as the bit name.  Each
            key is given a bit number ranging from 0..N-1.

        descr (:obj:`str`, :obj:`list`, optional):
            List of descriptions (or single discription) provided by
            :func:`info` for each bit.  No descriptions by default.

    Raises:
        ValueError:
            Raised if more than 64 bits are provided.
        TypeError:
            Raised if the provided `keys` do not have the correct type.

    Attributes:
        nbits (int):
            Number of bits
        bits (dict):
            A dictionary with the bit name and value
        descr (numpy.ndarray):
            List of bit descriptions
        max_value (int):
            The maximum valid bitmask value given the number of bits.
    """
    def __init__(self, keys, descr=None):

        _keys = keys if hasattr(keys, '__iter__') else [keys]
        _keys = numpy.atleast_1d(_keys).ravel()
        _descr = None if descr is None else numpy.atleast_1d(descr).ravel()

        if not numpy.all([isinstance(k, str) for k in _keys]):
            raise TypeError('Input keys must have string type.')
        if _descr is not None:
            if not all([isinstance(d, str) for d in _descr]):
                raise TypeError('Input keys must have string type.')
            if len(_descr) != len(_keys):
                raise ValueError('Number of listed descriptions not the same as number of keys.')

        # Do not allow for more that 64 bits
        if len(_keys) > 64:
            raise ValueError('Can only define up to 64 bits!')

        # Allow for multiple NULL keys; but check the rest for
        # uniqueness
        diff = set(_keys) - set(['NULL'])
        if len(diff) != numpy.unique(_keys[_keys != 'NULL']).size:
            raise ValueError('All input keys must be unique.')

        # Initialize the attributes
        self.nbits = len(_keys)
        self.bits = { k:i for i,k in enumerate(_keys) }
        self.max_value = (1 << self.nbits)-1
        self.descr = _descr
        
    def _prep_flags(self, flag):
        """Prep the flags for use."""
        # Flags must be a numpy array
        _flag = numpy.array(self.keys()) if flag is None else numpy.atleast_1d(flag).ravel()
        # NULL flags not allowed
        if numpy.any(_flag == 'NULL'):
            raise ValueError('Flag name NULL is not allowed.')
        # Flags should be among the bitmask keys
        if numpy.any([f not in self.keys() for f in _flag]):
            raise ValueError('Some bit names not recognized.')
        # Flags should be strings
        if numpy.any([ not isinstance(f, str) for f in _flag ]):
            raise TypeError('Provided bit names must be strings!')
        return _flag

    def keys(self):
        """
        Return a list of the bits; 'NULL' keywords are ignored.
        
        Returns:
            list: List of bit keywords.
        """
        return list(set(self.bits.keys())-set(['NULL']))

    def info(self):
        """
        Print the list of bits and, if available, their descriptions.
        """
        try:
            tr, tcols = numpy.array(os.popen('stty size', 'r').read().split()).astype(int)
            tcols -= int(tcols*0.1)
        except:
            tr = None
            tcols = None

        for k,v in sorted(self.bits.items(), key=lambda x:(x[1],x[0])):
            if k == 'NULL':
                continue
            print('         Bit: {0} = {1}'.format(k,v))
            if self.descr is not None:
                if tcols is not None:
                    print(textwrap.fill(' Description: {0}'.format(self.descr[v]), tcols))
                else:
                    print(' Description: {0}'.format(self.descr[v]))
            print(' ')

    def minimum_dtype(self, asuint=False):
        """
        Return the smallest int datatype that is needed to contain all
        the bits in the mask.  Output as an unsigned int if requested.

        Args:
            asuint (:obj:`bool`, optional):
                Return an unsigned integer type.  Signed types are
                returned by default.

        .. warning::
            uses int16 if the number of bits is less than 8 and
            asuint=False because of issue astropy.io.fits has writing
            int8 values.
        """
        if self.nbits < 8:
            return numpy.uint8 if asuint else numpy.int16
        if self.nbits < 16:
            return numpy.uint16 if asuint else numpy.int16
        if self.nbits < 32:
            return numpy.uint32 if asuint else numpy.int32
        return numpy.uint64 if asuint else numpy.int64

    def flagged(self, value, flag=None):
        """
        Determine if a bit is on in the provided bitmask value.  The
        function can be used to determine if any individual bit is on or
        any one of many bits is on.

        Args:
            value (int, array-like):
                Bitmask value.  It should be less than or equal to
                :attr:`max_value`; however, that is not checked.
            flag (str, array-like, optional):
                One or more bit names to check.  If None, then it checks
                if *any* bit is on.
        
        Returns:
            bool: Boolean flags that the provided flags (or any flag) is
            on for the provided bitmask value.  Shape is the same as
            `value`.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            TypeError: Raised if the provided *flag* does not contain
                one or more strings.
        """
        _flag = self._prep_flags(flag)

        out = value & (1 << self.bits[_flag[0]]) != 0
        if len(_flag) == 1:
            return out

        nn = len(_flag)
        for i in range(1,nn):
            out |= (value & (1 << self.bits[_flag[i]]) != 0)
        return out

    def flagged_bits(self, value):
        """
        Return the list of flagged bit names for a single bit value.

        Args:
            value (int):
                Bitmask value.  It should be less than or equal to
                :attr:`max_value`; however, that is not checked.
        
        Returns:
            list: List of flagged bit value keywords.

        Raises:
            KeyError:
                Raised by the dict data type if the input *flag* is not
                one of the valid :attr:`flags`.
            TypeError:
                Raised if the provided *flag* does not contain one or
                more strings.
        """
        if not numpy.issubdtype(type(value), numpy.integer):
            raise TypeError('Input must be a single integer.')
        if value <= 0:
            return []
        keys = numpy.array(self.keys())
        indx = numpy.array([1<<self.bits[k] & value != 0 for k in keys])
        return list(keys[indx])

    def toggle(self, value, flag):
        """
        Toggle a bit in the provided bitmask value.

        Args:
            value (int, array-like):
                Bitmask value.  It should be less than or equal to
                :attr:`max_value`; however, that is not checked.
            flag (str, array-like):
                Bit name(s) to toggle.

        Returns:
            array-like: New bitmask value after toggling the selected
            bit.

        Raises:
            ValueError
            KeyError:

                Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.

            Exception: Raised if the provided *flag* is not a string.
        """ 
        if flag is None:
            raise ValueError('Provided bit name cannot be None.')

        _flag = self._prep_flags(flag)

        out = value ^ (1 << self.bits[_flag[0]])
        if len(_flag) == 1:
            return out

        nn = len(_flag)
        for i in range(1,nn):
            out ^= (1 << self.bits[_flag[i]])
        return out

    def turn_on(self, value, flag):
        """
        Ensure that a bit is turned on in the provided bitmask value.

        Args:
            value (uint or array): Bitmask value.  It should be less
                than or equal to :attr:`max_value`; however, that is not
                checked.
            flag (list, numpy.ndarray, or str): Bit name(s) to turn on.
        
        Returns:
            uint: New bitmask value after turning on the selected bit.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            Exception: Raised if the provided *flag* is not a string.
        """
        if flag is None:
            raise ValueError('Provided bit name cannot be None.')

        _flag = self._prep_flags(flag)

        out = value | (1 << self.bits[_flag[0]])
        if len(_flag) == 1:
            return out

        nn = len(_flag)
        for i in range(1,nn):
            out |= (1 << self.bits[_flag[i]])
        return out

    def turn_off(self, value, flag):
        """
        Ensure that a bit is turned off in the provided bitmask value.

        Args:
            value (uint or array): Bitmask value.  It should be less
                than or equal to :attr:`max_value`; however, that is not
                checked.
            flag (list, numpy.ndarray, or str): Bit name(s) to turn off.
        
        Returns:
            uint: New bitmask value after turning off the selected bit.

        Raises:
            KeyError: Raised by the dict data type if the input *flag*
                is not one of the valid :attr:`flags`.
            Exception: Raised if the provided *flag* is not a string.
        """
        if flag is None:
            raise ValueError('Provided bit name cannot be None.')

        _flag = self._prep_flags(flag)

        out = value & ~(1 << self.bits[_flag[0]])
        if len(_flag) == 1:
            return out

        nn = len(_flag)
        for i in range(1,nn):
            out &= ~(1 << self.bits[_flag[i]])
        return out

    def consolidate(self, value, flag_set, consolidated_flag):
        """
        Consolidate a set of flags into a single flag.
        """
        indx = self.flagged(value, flag=flag_set)
        value[indx] = self.turn_on(value[indx], consolidated_flag)
        return value

