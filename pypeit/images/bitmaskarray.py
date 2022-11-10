"""
General utility for bit-based masking.
"""
from IPython import embed

import numpy as np

from pypeit.datamodel import DataContainer


class BitMaskArray(DataContainer):
    version = None

    datamodel = {'mask': dict(otype=np.ndarray, atype=np.integer, descr='Bitmask values')}

    bitmask = None

    def __init__(self, shape, asuint=False):
        # Instantiate as an empty DataContainer
        super().__init__()

        self.lower_keys = [k.lower() for k in self.bit_keys()]
        self.mask = np.zeros(shape, dtype=self.bitmask.minimum_dtype(asuint=asuint))

    def __getattr__(self, item):
        try:
            i = self.lower_keys.index(item.lower())
        except ValueError:
            return super().__getattr__(item)
        # TODO: This creates a new numpy array every time it is called, which
        # could be slow for large arrays.  Instead, we might want to do some
        # clever lazy-loading to make this faster.
        return self.flagged(flag=list(self.bit_keys())[i])

    def _init_internals(self):
        """
        Initialize attributes that are not part of the datamodel.
        """
        self.lower_keys = None

    # TODO: This loses the description of the bits.  Might be better to override
    # to_hdu; although this gets sticky trying to figure out which hdu has the
    # mask...  Leaving this way for now.
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

    # NOTE: This function cannot be called keys because that would override the
    # DataContainer base-class function!
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

    def flagged(self, flag=None):
        """
        Determine if a bit is on in the provided bitmask value.  The
        function can be used to determine if any individual bit is on or
        any one of many bits is on.

        Args:
            flag (str, array-like, optional):
                One or more bit names to check.  If None, then it checks
                if *any* bit is on.
        
        Returns:
            bool: Boolean flags that the provided flags (or any flag) is
            on for the provided bitmask value.  Shape is the same as
            `value`.
        """
        return self.bitmask.flagged(self.mask, flag=flag)

    def flagged_bits(self, index):
        """
        Return the list of flagged bit names for a single bit value.

        Args:
            index (:obj:`tuple`):
                Tuple with the indices in the array.
        
        Returns:
            list: List of flagged bit value keywords.
        """
        return self.bitmask.flagged_bits(self.mask[index])

    def toggle(self, select, flag):
        """
        Toggle bits for selected array elements.

        Args:
            select (:obj:`tuple`, :obj:`slice`, `numpy.ndarray`_):
                Object used to select elements of the mask array to at which to
                toggle the provided bit flags.  I.e., for the internal
                :attr:`mask`, ``mask[select]`` must be a valid (fancy indexing)
                operation.
            flag (:obj:`str`, array-like):
                Bit name(s) to toggle.
        """
        self.mask[select] = self.bitmask.toggle(self.mask[select], flag)

    def turn_on(self, select, flag):
        """
        Ensure that a bit is turned on for the selected elements.

        Args:
            select (:obj:`tuple`, :obj:`slice`, `numpy.ndarray`_):
                Object used to select elements of the mask array to at which to
                turn on the provided bit flags.  I.e., for the internal
                :attr:`mask`, ``mask[select]`` must be a valid (fancy indexing)
                operation.
            flag (:obj:`str`, array-like):
                Bit name(s) to turn on.
        """
        self.mask[select] = self.bitmask.turn_on(self.mask[select], flag)

    def turn_off(self, select, flag):
        """
        Ensure that a bit is turned off in the provided bitmask value.

        Args:
            select (:obj:`tuple`, :obj:`slice`, `numpy.ndarray`_):
                Object used to select elements of the mask array to at which to
                turn off the provided bit flags.  I.e., for the internal
                :attr:`mask`, ``mask[select]`` must be a valid (fancy indexing)
                operation.
            flag (str, array-like):
                Bit name(s) to turn off.
        """
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


