
Bitmasks allow you to define a set of bit values signified by strings,
and then toggle and interpret bits held by a `numpy.ndarray`_.  For
example, say you're processing an image and you want to set up a set of
bits that indicate that the pixel is part of a bad-pixel mask, has a
cosmic ray, or is saturated.  You can define the following::

    from pypeit.bitmask import BitMask

    bits = {'BPM':'Pixel is part of a bad-pixel mask',
            'COSMIC':'Pixel is contaminated by a cosmic ray',
            'SATURATED':'Pixel is saturated.'}
    image_bm = BitMask(list(bits.keys()), descr=list(bits.values()))

.. note::

    Consistency in the order of the dictionary keywords is *critical*
    to the repeatability of the :class:`~pypeit.bitmask.BitMask`
    instances. The above is possible because :obj:`dict` objects
    automatically maintain the order of the provided keywords since
    Python 3.7, the minimum required version for ``PypeIt``.
    
Or, better yet, define a derived class::

    from pypeit.bitmask import BitMask

    class ImageBitMask(BitMask):
        def __init__(self):
            bits = {'BPM':'Pixel is part of a bad-pixel mask',
                    'COSMIC':'Pixel is contaminated by a cosmic ray',
                    'SATURATED':'Pixel is saturated.'}
            super(ImageBitMask, self).__init__(list(bits.keys()), descr=list(bits.values()))

    image_bm = ImageBitMask()

In either case, you can see the list of bits and their bit numbers by
running::

    >>> image_bm.info()
             Bit: BPM = 0
     Description: Pixel is part of a bad-pixel mask

             Bit: COSMIC = 1
     Description: Pixel is contaminated by a cosmic ray

             Bit: SATURATED = 2
     Description: Pixel is saturated.
    >>> image_bm.bits
    {'BPM': 0, 'COSMIC': 1, 'SATURATED': 2}
    >>> image_bm.keys()
    ['BPM', 'COSMIC', 'SATURATED']

Now you can define a `numpy.ndarray`_ to hold the mask value for each
image pixel; the :func:`~pypeit.bitmask.BitMask.minimum_dtype`
returns the the smallest data type required to represent the list of
defined bits. The maximum number of bits that can be defined is 64.
Assuming you have an image ``img``::

    import numpy
    mask = numpy.zeros(img.shape, dtype=image_bm.minimum_dtype())

Assuming you have boolean or integer arrays that identify pixels to
mask, you can turn on the mask bits as follows::

    mask[cosmics_indx] = image_bm.turn_on(mask[cosmics_indx], 'COSMIC')
    mask[saturated_indx] = image_bm.turn_on(mask[saturated_indx], 'SATURATED')

or make sure certain bits are off::

    mask[not_a_cosmic] = image_bm.turn_off(mask[not_a_cosmic], 'COSMIC')

The form of these methods is such that the array passed to the method
are not altered.  Instead the altered bits are returned, which is why
the lines above have the form ``m = bm.turn_on(m, flag)``.

Some other short usage examples:

    - To find which flags are set for a single value::
        
        image_bm.flagged_bits(mask[0,10])

    - To find the list of unique flags set for any pixel::

        unique_flags = numpy.sort(numpy.unique(numpy.concatenate(
                            [image_bm.flagged_bits(b) for b in numpy.unique(mask)]))).tolist()

    - To get a boolean array that selects pixels with one or more
      mask bits::

        cosmics_indx = image_bm.flagged(mask, flag='COSMIC')
        all_but_bpm_indx = image_bm.flagged(mask, flag=['COSMIC', 'SATURATED'])
        any_flagged = image_bm.flagged(mask)

    - To construct masked arrays, following from the examples above::

        masked_img = numpy.ma.MaskedArray(img, mask=image_bm.flagged(mask))

:class:`~pypeit.bitmask.BitMask` objects can be defined
programmatically, as shown above for the ``ImageBitMask`` derived class,
but they can also be defined by reading formatted files.  The current
options are:

    #. Fits headers: There are both reading and writing methods for
       bitmask I/O using `astropy.io.fits.Header`_ objects.  Using the
       ``ImageBitMask`` class as an example::
       
            >>> from astropy.io import fits
            >>> hdr = fits.Header()
            >>> image_bm = ImageBitMask()
            >>> image_bm.to_header(hdr)
            >>> hdr
            BIT0    = 'BPM     '           / Pixel is part of a bad-pixel mask
            BIT1    = 'COSMIC  '           / Pixel is contaminated by a cosmic ray
            BIT2    = 'SATURATED'          / Pixel is saturated.
            >>> copy_bm = BitMask.from_header(hdr)

