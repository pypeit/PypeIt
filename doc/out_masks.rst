
:orphan:

.. include:: include/links.rst

.. _out_masks:

===============
Output Bitmasks
===============

PypeIt uses bitmasks to flag data for various reasons. For a primer on the
general use of bitmasks, see the `SDSS bitmasks`_ primer.

For the underlying bitmask utility class used by PypeIt, see :doc:`bitmasks`.

----

Spec2D Bits
===========

In terms of the output data, the main bitmask you will encounter is
in the :ref:`spec-2d-output`. The ``BPMMASK`` extension in these files
contain the bad-pixel bitmask values accrued during the processing
and reduction of the data. The main object that manages the bitmasks
is :class:`~pypeit.images.imagebitmask.ImageBitMask`, which defines
the following bits:

.. include:: include/imagebitmask_table.rst

Interpreting the bit values
---------------------------

You can construct the appropriate bitmask in two ways:

    - (Recommended) Instantiate the bitmask directly:

      .. code-block:: python

        from pypeit.images.imagebitmask import ImageBitMask
        bm = ImageBitMask()

    - Specifically for the ``spec2d*`` files, the names of the bits
      and their order is saved to the header. You can instantiate a
      generic :class:`~pypeit.bitmask.BitMask` from the header;
      however, it's not clear how long this behavior will remain. To
      instantiate the relevant :class:`~pypeit.bitmask.BitMask` from
      the header:

      .. code-block:: python

        from astropy.io import fits
        from pypeit.bitmask import BitMask
        hdu = fits.open('spec2d_DE.20100913.22358-CFHQS1_DEIMOS_20100913T061231.334.fits')
        bm = BitMask(hdu[1].header['IMGBITM'].split(','))


With the :class:`~pypeit.bitmask.BitMask` or
:class:`~pypeit.images.imagebitmask.ImageBitMask` instantiated, you
can interpret the meaning of the ``BPMMASK`` extensions as follows:

    - Use the :func:`~pypeit.bitmask.BitMask.flagged` method to
      produce a boolean array that selects pixels flagged with a
      given bit:

      .. code-block:: python

        extract_flag = bm.flagged(hdu['DET01-BPMMASK'].data, flag='EXTRACT')

    - Use the same method to select multiple bits:

      .. code-block:: python

        process_flag = bm.flagged(hdu['DET01-BPMMASK'].data, flag=['BPM', 'CR', 'SATURATION'])

    - Or select bits that are flagged for *any* reason:

      .. code-block:: python

        # You can use the `flagged` method for this
        any_flag = bm.flagged(hdu['DET01-BPMMASK'].data)
        # or equivalently
        _any_flag = hdu['DET01-BPMMASK'].data > 0

    - To get the human-readable reason that any given value is flagged, use the
      :func:`~pypeit.bitmask.BitMask.flagged_bits` method.  Currently this can
      only be used with a *single* bit value:

      .. code-block:: python

        print(bm.flagged_bits(hdu['DET01-BPMMASK'].data[0,0])) 

