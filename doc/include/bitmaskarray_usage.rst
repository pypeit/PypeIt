
:class:`~pypeit.images.bitmaskarray.BitMaskArray` objects combine the `numpy.ndarray`_
that holds the bit values and the :class:`~pypeit.bitmask.BitMask` object used
to interpret them into a new subclass of
:class:`~pypeit.datamodel.DataContainer`.  The following builds on the example
uses of :class:`~pypeit.bitmask.BitMask` objects (see :ref:`bitmasks`).

Defining a new subclass
+++++++++++++++++++++++

Given the definition of ``ImageBitMask``, we can implement a relevant
:class:`~pypeit.images.bitmaskarray.BitMaskArray` subclass as follows:

.. code-block:: python

    from pypeit.images.bitmaskarray import BitMaskArray

    class ImageBitMaskArray(BitMaskArray):
        version = '1.0'
        bitmask = ImageBitMask()

The remaining functionality below is all handled by the base class.

To instantiate a new 2D mask array that is 5 pixels on a side:

.. code-block:: python

    shape = (5,5)
    mask = ImageBitMaskArray(shape)

You can access the bit flag names using:

.. code-block:: python

    >>> mask.bit_keys()
    ['BPM', 'COSMIC', 'SATURATED']

Bit access
++++++++++

You can flag bits using :func:`~pypeit.images.bitmaskarray.BitMaskArray.turn_on`.  For
example, the following code flags the center column of the image as being
part of the detector bad-pixel mask:

.. code-block:: python

    import numpy as np
    mask.turn_on('BPM', select=np.s_[:,2])

The ``select`` argument to
:func:`~pypeit.images.bitmaskarray.BitMaskArray.turn_on` can be *anything* that
is appropriately interpreted as slicing a `numpy.ndarray`_.  That is,
``arr[select]`` must be valid, where ``arr`` is the internal array held by
``mask``.

Similarly, you can flag a pixel with a cosmic ray:

.. code-block:: python

    mask.turn_on('COSMIC', select=(0,0))

or multiple pixels as being saturated:

.. code-block:: python

    mask.turn_on('SATURATED', select=([0,1,-1,-1],[0,0,-1,-2]))

and you can simultaneously flag pixels for multiple reasons:

.. code-block:: python

    mask.turn_on(['COSMIC', 'SATURATED'], select=([-1,-1],[0,1]))


The mask values themselves are accessed using the ``mask`` attribute:

.. code-block:: python
    
    >>> mask.mask
    array([[6, 0, 1, 0, 0],
           [4, 0, 1, 0, 0],
           [0, 0, 1, 0, 0],
           [0, 0, 1, 0, 0],
           [6, 6, 1, 4, 4]], dtype=int16)

However, more usefully, you can obtain boolean arrays that select pixels flagged
by one or more flags:

.. code-block:: python

    >>> mask.flagged(flag='SATURATED')
    array([[ True, False, False, False, False],
           [ True, False, False, False, False],
           [False, False, False, False, False],
           [False, False, False, False, False],
           [ True,  True, False,  True,  True]])

    >>> mask.flagged(flag=['BPM', 'SATURATED'])
    array([[ True, False,  True, False, False],
           [ True, False,  True, False, False],
           [False, False,  True, False, False],
           [False, False,  True, False, False],
           [ True,  True,  True,  True,  True]])

If you want to select all pixels that are **not** flagged by a given flag, you
can use the ``invert`` option in 
:func:`~pypeit.images.bitmaskarray.BitMaskArray.flagged`:

.. code-block:: python

    >>> gpm = mask.flagged(flag='BPM', invert=True)
    >>> gpm
    array([[ True,  True, False,  True,  True],
           [ True,  True, False,  True,  True],
           [ True,  True, False,  True,  True],
           [ True,  True, False,  True,  True],
           [ True,  True, False,  True,  True]])

For individual flags, there is also convenience functionality that allows you to
access a boolean array as if it were an attribute of the object:

.. code-block:: python

    >>> mask.bpm
    array([[False, False,  True, False, False],
           [False, False,  True, False, False],
           [False, False,  True, False, False],
           [False, False,  True, False, False],
           [False, False,  True, False, False]])
    >>> mask.saturated
    array([[ True, False, False, False, False],
           [ True, False, False, False, False],
           [False, False, False, False, False],
           [False, False, False, False, False],
           [ True,  True, False,  True,  True]])

This convenience operation is identical to calling
:func:`~pypeit.images.bitmaskarray.BitMaskArray.flagged` for the indicated bit.
However ``bpm`` is **not** an array that can be used to change the value of the
bits themselves:

.. code-block:: python

    >>> indx = np.zeros(shape, dtype=bool)
    >>> indx[2,3] = True
    >>> mask.bpm = indx # Throws an AttributeError

Instead, you must use the bit toggling functions provided by the class:
:func:`~pypeit.images.bitmaskarray.BitMaskArray.turn_on`,
:func:`~pypeit.images.bitmaskarray.BitMaskArray.turn_off`, or
:func:`~pypeit.images.bitmaskarray.BitMaskArray.toggle`.

.. tip::

    Every time :func:`~pypeit.images.bitmaskarray.BitMaskArray.flagged` is
    called, a new array is created.  If you need to access to the result of the
    function multiple times without changing the flags, you're better of
    assigning the result to a new array and then using that array so that you're
    not continually allocating and deallocating memory (even within the context
    of how this is done within python).

Input/Output
++++++++++++

As a subclass of :class:`~pypeit.datamodel.DataContainer`, you can save and read
the bitmask data to/from files:

.. code-block:: python

    >>> mask.to_file('mask.fits')
    >>> _mask = ImageBitMaskArray.from_file('mask.fits')
    >>> np.array_equal(mask.mask, _mask.mask)
    True

In addition to the mask data, the bit flags and values are also written to the
header; see the ``BIT*`` entries in the header below:

.. code-block:: python

    >>> from astropy.io import fits
    >>> hdu = fits.open('mask.fits')
    >>> hdu.info()
    Filename: mask.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU      13   ()
      1  MASK          1 ImageHDU        22   (5, 5)   int16
    >>> hdu['MASK'].header
    XTENSION= 'IMAGE   '           / Image extension
    BITPIX  =                   16 / array data type
    NAXIS   =                    2 / number of array dimensions
    NAXIS1  =                    5
    NAXIS2  =                    5
    PCOUNT  =                    0 / number of parameters
    GCOUNT  =                    1 / number of groups
    VERSPYT = '3.9.13  '           / Python version
    VERSNPY = '1.22.3  '           / Numpy version
    VERSSCI = '1.8.0   '           / Scipy version
    VERSAST = '5.0.4   '           / Astropy version
    VERSSKL = '1.0.2   '           / Scikit-learn version
    VERSPYP = '1.10.1.dev260+g32de3d6d4' / PypeIt version
    DATE    = '2022-11-10'         / UTC date created
    DMODCLS = 'ImageBitMaskArray'  / Datamodel class
    DMODVER = '1.0     '           / Datamodel version
    BIT0    = 'BPM     '
    BIT1    = 'COSMIC  '
    BIT2    = 'SATURATED'
    EXTNAME = 'MASK    '           / extension name
    CHECKSUM= 'APGODMFOAMFOAMFO'   / HDU checksum updated 2022-11-10T13:10:27
    DATASUM = '1245200 '           / data unit checksum updated 2022-11-10T13:10:27

.. note::

    Currently, when loading a mask, the bit names in the header of the output
    file are **not** checked against the bitmask definition in the code itself.
    This kind of version control should be handled using the ``version``
    attribute of the class.  I.e., anytime the flags in the bitmask are changed,
    the developers should bump the class version.

