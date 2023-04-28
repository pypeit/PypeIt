"""
General utility for bit-based masking.

Class usage examples
====================

.. include:: ../include/bitmaskarray_usage.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from IPython import embed

import numpy as np

from pypeit.datamodel import DataContainer
from pypeit import msgs


class BitMaskArray(DataContainer):
    """
    Object that holds both the mask data and the mask interpretation object.

    This is an abstract class that should not be directly instantiated.

    Args:
        shape (:obj:`tuple`):
            Shape of the mask to create.
        asuint (:obj:`bool`, optional):
            When setting the data-type for the mask array (see
            :func:`~pypeit.bitmask.BitMask.minimum_dtype`), use an *unsigned*
            integer instead of a signed integer (e.g., ``uint16`` instead of
            ``int16``).
    """

    version = None
    """
    DataContainer version.  Must be defined by the subclass.
    """

    datamodel = {'mask': dict(otype=np.ndarray, atype=np.integer, descr='Bitmask values')}
    """
    Datamodel is simple, containing only the mask array.
    """

    internals = ['lower_keys']

    bitmask = None
    """
    :class:`~pypeit.bitmask.BitMask` object used to interpret the bit array.
    Must be defined by the subclass.  When defining subclasses, note that the
    bitmask flags *must* be case-insensitive strings.
    """

    def __init__(self, shape, asuint=False):
        # Instantiate as an empty DataContainer
        super().__init__()
        self._set_keys()
        self.mask = np.zeros(shape, dtype=self.bitmask.minimum_dtype(asuint=asuint))

    def _set_keys(self):
        """
        Set :attr:`lower_keys`, which are needed for the bit access convenience
        method.
        """
        # Check the bitmask
        keys = self.bit_keys()
        if any([not isinstance(k, str) for k in keys]):
            msgs.error(f'CODING ERROR: {self.bitmask.__class__.__name__} must only contain '
                       'string bit flags.')

        self.lower_keys = [k.lower() for k in keys]
        if len(np.unique(self.lower_keys)) != len(keys):
            msgs.error('CODING ERROR: All bitmask keys must be case-insensitive and unique: '
                       f'{keys}')

    def __getattr__(self, item):
        """
        Override the attribute access to allow for immediate construction and
        access to boolean arrays that select array elements with the provided flag.

        For example, if ``'BPM'`` is a flag in the :attr:`bitmask`, and ``mask``
        is an instance of the the relevant subclass, ``mask.bpm`` is identical
        to calling ``mask.flagged(flag='BPM')``.

        .. warning::

            Using this functionality creates a new numpy array *every time it is
            called*.  This can be slow for large arrays.  This means that, if
            you need to access the boolean array for a given flag multiple
            times, you should set it to a new object (e.g., ``maskbpm =
            mask.bpm``)!

        If the attribute (``item``) requested is *not* one of the bitmask flags,
        the :class:`~pypeit.datamodel.DataContainer` base class function is
        called.

        Args:
            item (object):
                The attribute being accessed.
        """
        try:
            i = self.lower_keys.index(item.lower())
        except ValueError:
            return super().__getattr__(item)
        # TODO: This creates a new numpy array every time it is called, which
        # could be slow for large arrays.  Instead, we might want to do some
        # clever lazy-loading to make this faster.
        return self.flagged(flag=list(self.bit_keys())[i])

    def __getitem__(self, item):
        """Allow direct access to the mask."""
        try:
            return self.mask[item]
        except:
            return super().__getitem__(item)

    def __setitem__(self, item, value):
        try:
            self.mask[item] = value
        except:
            super().__setitem__(item, value)

    def __or__(self, other):
        """Override or operation for mask."""
        _self = self.copy()
        _self.mask |= other.mask
        return _self

    def __and__(self, other):
        """Override and operation for mask."""
        _self = self.copy()
        _self.mask &= other.mask
        return _self

    # TODO: This loses the description of the bits.  Might be better to override
    # to_hdu; although this gets sticky trying to figure out which hdu has the
    # mask...  Leaving it this way for now.
    def _bundle(self):
        """
        Override the base-class bundle method so that the bitmask keys can be
        added to the header.

        Returns:
            :obj:`list`: List of dictionaries indicating what should be written
            to a file.
        """
        d = super()._bundle()
        d[0].update(self.bitmask.to_dict())
        return d

    def copy(self):
        """Create a deep copy."""
        _self = super().__new__(self.__class__)
        DataContainer.__init__(_self)
        _self._set_keys()
        _self.mask = self.mask.copy()
        return _self

    # NOTE: This function cannot be called "keys" because that would override
    # the DataContainer base-class function!
    def bit_keys(self):
        """
        Return a list of the bit keywords.
        """
        return self.bitmask.keys()

    @property
    def bits(self):
        """Return the bit dictionary."""
        return self.bitmask.bits

    @property
    def shape(self):
        """Return the shape of the internal array."""
        return self.mask.shape

    def info(self):
        """
        Print the list of bits and, if available, their descriptions.
        """
        self.bitmask.info()

    def flagged(self, flag=None, invert=False):
        """
        Determine if a bit is on in the provided bitmask value.  The
        function can be used to determine if any individual bit is on or
        any one of many bits is on.

        Args:
            flag (:obj:`str`, array-like, optional):
                One or more bit names to check.  If None, then it checks
                if *any* bit is on.
            invert (:obj:`bool`, optional):
                Invert the boolean such that unflagged pixels are True and
                flagged pixels are False.
        
        Returns:
            `numpy.ndarray`_:  Boolean array indicating where the internal
            bitmask is flagged by the selected bits.  Flagged values are True,
            unflagged values are False.  If ``invert`` is True, this is
            reversed.
        """
        indx = self.bitmask.flagged(self.mask, flag=flag)
        return np.logical_not(indx) if invert else indx

    def flagged_bits(self, index):
        """
        Return the list of flagged bit names for a single bit value.

        Args:
            index (:obj:`tuple`):
                Tuple with the indices in the array.
        
        Returns:
            :obj:`list`: List of flagged bit value keywords.
        """
        return self.bitmask.flagged_bits(self.mask[index])

    def toggle(self, flag, select=None):
        """
        Toggle bits for selected array elements.

        Args:
            flag (:obj:`str`, array-like):
                Bit name(s) to toggle.
            select (:obj:`tuple`, :obj:`slice`, `numpy.ndarray`_, optional):
                Object used to select elements of the mask array to at which to
                toggle the provided bit flags.  I.e., for the internal
                :attr:`mask`, ``mask[select]`` must be a valid (fancy indexing)
                operation.  If None, the bit is toggled for the full mask!
        """
        if select is None:
            select = np.s_[...]
        self.mask[select] = self.bitmask.toggle(self.mask[select], flag)

    def turn_on(self, flag, select=None):
        """
        Ensure that a bit is turned on for the selected elements.

        Args:
            flag (:obj:`str`, array-like):
                Bit name(s) to turn on.
            select (:obj:`tuple`, :obj:`slice`, `numpy.ndarray`_, optional):
                Object used to select elements of the mask array to at which to
                turn on the provided bit flags.  I.e., for the internal
                :attr:`mask`, ``mask[select]`` must be a valid (fancy indexing)
                operation.  If None, the bit is turned on for the full mask!
        """
        if select is None:
            select = np.s_[...]
        self.mask[select] = self.bitmask.turn_on(self.mask[select], flag)

    def turn_off(self, flag, select=None):
        """
        Ensure that a bit is turned off in the provided bitmask value.

        Args:
            flag (:obj:`str`, array-like):
                Bit name(s) to turn off.
            select (:obj:`tuple`, :obj:`slice`, `numpy.ndarray`_, optional):
                Object used to select elements of the mask array to at which to
                turn off the provided bit flags.  I.e., for the internal
                :attr:`mask`, ``mask[select]`` must be a valid (fancy indexing)
                operation.  If None, the bit is turned off for the full mask!
        """
        if select is None:
            select = np.s_[...]
        self.mask[select] = self.bitmask.turn_off(self.mask[select], flag)

    def consolidate(self, flag_set, consolidated_flag):
        """
        Consolidate a set of flags into a single flag.

        That is, any bit flagged with any of the flags provided by ``flag_set``
        will also be flagged by ``consolidate_flag`` after executing this
        function.

        Args:
            flag_set (:obj:`str`, array-like):
                List of flags that are consolidated into a single flag.
            consolidated_flag (:obj:`str`):
                Consolidated flag name.
        """
        self.turn_on(self.flagged(flag=flag_set), consolidated_flag)


