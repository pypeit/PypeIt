
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

Viewing and interpreting the bit mask image (python)
----------------------------------------------------

To access the bitmask array, we recommend using
:class:`~pypeit.spec2dobj.Spec2DObj` directly:

.. code-block:: python

    from pathlib import Path
    from pypeit.spec2dobj import Spec2DObj

    f = Path('spec2d_b27-J1217p3905_KASTb_20150520T045733.560.fits').resolve()
    detname = 'DET01'
    spec2d = Spec2DObj.from_file(f, detname)
    mask = spec2d.bpmmask

The bitmask array is held by the ``spec2d.bpmmask`` attribute.

But you can read it directly from the spec2d file:

.. code-block:: python

    from astropy.io import fits
    from pypeit.images.imagebitmask import ImageBitMaskArray

    f = Path('spec2d_b27-J1217p3905_KASTb_20150520T045733.560.fits').resolve()
    detname = 'DET01'
    hdu = fits.open(f)
    mask = ImageBitMaskArray.from_array(hdu[f'{detname}-BPMMASK'].data)

To show all the bit values directly:

.. code-block:: python

    from matplotlib import pyplot as plt

    plt.imshow(mask, origin='lower', interpolation='nearest')
    plt.show()

However, this isn't necessarily as useful as creating boolean arrays that
identify which pixels are flagged due to one or more reasons.

E.g., to show all pixels flagged for having cosmic rays:

.. code-block:: python

    plt.imshow(mask.flagged(flag='CR').astype(int), origin='lower', interpolation='nearest')
    plt.show()

or as being part of the instrument-specific bad-pixel mask *or* not part of any slit:

.. code-block:: python

    plt.imshow(mask.flagged(flag=['BPM', 'OFFSLITS']).astype(int),
               origin='lower', interpolation='nearest')
    plt.show()

You can also use the :ref:`pypeit_show_2dspec` script to include an image that
shows the full mask or an image that selects specific flags.

To print the human-readable reason(s) any given value is flagged:

.. code-block:: python

    coo = (0,0)   # Tuple with the 2D coordinate of the pixel
    print(mask.flagged_bits(coo))

.. _pypeit_print_bpm:

pypeit_print_bpm
----------------

This simple executable allows you to effectively do the above via a command-line
script.  The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_print_bpm.rst

A typical call and its output looks like this:

.. code-block:: bash

    % pypeit_print_bpm 23
    [INFO]    :: Using the default PypeIt bad pixel mask.
    [INFO]    :: The bad pixel mask value (23) corresponds to the following:

                 * BPM        : Component of the instrument-specific bad pixel mask
                 * CR         : Cosmic ray detected
                 * SATURATION : Saturated pixel
                 * OFFSLITS   : Pixel does not belong to any slit

    [INFO]    :: Please see the following website for more information:
                 https://pypeit.readthedocs.io/en/release/out_masks.html


