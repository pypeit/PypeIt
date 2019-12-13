# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides a set of I/O routines.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import sys
import warnings
import gzip
import shutil
from packaging import version

import numpy

from astropy.io import fits

# These imports are largely just to make the versions available for
# writing to the header. See `initialize_header`
import scipy
import astropy
import sklearn
import pypeit
import time
from IPython import embed

def init_record_array(shape, dtype):
    r"""
    Utility function that initializes a record array using a provided
    input data type.  For example::

        dtype = [ ('INDX', numpy.int, (2,) ),
                  ('VALUE', numpy.float) ]

    Defines two columns, one named `INDEX` with two integers per row and
    the one named `VALUE` with a single float element per row.  See
    `numpy.recarray`_.
    
    Args:
        shape (:obj:`int`, :obj:`tuple`):
            Shape of the output array.
        dtype (:obj:`list`):
            List of the tuples that define each element in the record
            array.

    Returns:
        `numpy.recarray`: Zeroed record array
    """
    return numpy.zeros(shape, dtype=dtype).view(numpy.recarray)


def rec_to_fits_type(rec_element):
    """
    Return the string representation of a fits binary table data type
    based on the provided record array element.
    """
    n = 1 if len(rec_element[0].shape) == 0 else rec_element[0].size
    if rec_element.dtype == numpy.bool:
        return '{0}L'.format(n)
    if rec_element.dtype == numpy.uint8:
        return '{0}B'.format(n)
    if rec_element.dtype == numpy.int16 or rec_element.dtype == numpy.uint16:
        return '{0}I'.format(n)
    if rec_element.dtype == numpy.int32 or rec_element.dtype == numpy.uint32:
        return '{0}J'.format(n)
    if rec_element.dtype == numpy.int64 or rec_element.dtype == numpy.uint64:
        return '{0}K'.format(n)
    if rec_element.dtype == numpy.float32:
        return '{0}E'.format(n)
    if rec_element.dtype == numpy.float64:
        return '{0}D'.format(n)
    
    # If it makes it here, assume its a string
    l = int(rec_element.dtype.str[rec_element.dtype.str.find('U')+1:])
#    return '{0}A'.format(l) if n==1 else '{0}A{1}'.format(l*n,l)
    return '{0}A'.format(l*n)


def rec_to_fits_col_dim(rec_element):
    """
    Return the string representation of the dimensions for the fits
    table column based on the provided record array element.

    The shape is inverted because the first element is supposed to be
    the most rapidly varying; i.e. the shape is supposed to be written
    as row-major, as opposed to the native column-major order in python.
    """
    return None if len(rec_element[0].shape) == 1 else str(rec_element[0].shape[::-1])


def rec_to_bintable(arr, name=None, hdr=None):
    """
    Construct an `astropy.io.fits.BinTableHDU` from a record array.

    Args:
        arr (`numpy.recarray`):
            The data array to write to a binary table.
        name (:obj:`str`, optional):
            The name for the binary table extension.
        hdr (`astropy.io.fits.Header`, optional):
            Header for the BinTableHDU extension.
    
    Returns:
        `astropy.io.fits.BinTableHDU`: The binary fits table that can be
        included in an `astropy.io.fits.HDUList` and written to disk.
    """
    return fits.BinTableHDU.from_columns([fits.Column(name=n,
                                                      format=rec_to_fits_type(arr[n]),
                                                      dim=rec_to_fits_col_dim(arr[n]),
                                                      array=arr[n])
                                            for n in arr.dtype.names], name=name, header=hdr)


def compress_file(ifile, overwrite=False, rm_original=True):
    """
    Compress a file using gzip package.
    
    Args:
        ifile (:obj:`str`):
            Name of file to compress.  Output file with have the same
            name with '.gz' appended.
        overwrite (:obj:`bool`, optional):
            Overwrite any existing file.
        rm_original (:obj:`bool`, optional):
            The method writes the compressed file such that both the
            uncompressed and compressed file will exist when the
            compression is finished.  If this is True, the original
            (uncompressed) file is removed.

    Raises:
        ValueError:
            Raised if the file name already has a '.gz' extension or if
            the file exists and `overwrite` is False.
    """
    # Nominally check if the file is already compressed
    if ifile.split('.')[-1] == 'gz':
        raise ValueError('File appears to already have been compressed! {0}'.format(ifile))

    # Construct the output file name and check if it exists
    ofile = '{0}.gz'.format(ifile)
    if os.path.isfile(ofile) and not overwrite:
        raise FileExistsError('{0} exists! To overwrite, set overwrite=True.'.format(ofile))

    # Compress the file
    with open(ifile, 'rb') as f_in:
        with gzip.open(ofile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    if rm_original:
        # Remove the uncompressed file
        os.remove(ifile)


def parse_hdr_key_group(hdr, prefix='F'):
    """
    Parse a group of fits header values grouped by a keyword prefix.

    If the prefix is 'F', the header keywords are expected to be, e.g.,
    'F1', 'F2', 'F3', etc.  The list of values returned are then, e.g.,
    [ hdr['F1'], hdr['F2'], hdr['F3'], ... ].  The function performs no
    retyping, so the values in the returned list have whatever type they
    had in the fits header.

    Args:
        hdr (`fits.Header`):
            Astropy Header object
        prefix (:obj:`str`, optional):
            The prefix used for the header keywords.
        
    Returns:
        list: The list of header values ordered by their keyword index.
    """
    values = {}
    for k, v in hdr.items():
        # Check if this header keyword starts with the required prefix
        if k[:len(prefix)] == prefix:
            try:
                # Try to convert the keyword without the prefix into an
                # integer; offseting to an array index.
                i = int(k[len(prefix):])-1
            except ValueError:
                # Assume the value is some other random keyword that
                # starts with the prefix but isn't part of the keyword
                # group
                continue
            # Assume we've found a valid value; assign it and keep
            # trolling the header
            values[i] = v

    # Convert from dictionary with integer keys to an appropriately
    # sorted list
    # JFH try/except added to deal with cases where no header cards were found with the desired
    # prefix. I'm not sure returning an empty list is the desired behavior, but None will cause
    # downstream things to crash.
    try:
        return [values[i] for i in range(max(values.keys())+1)]
    except:
        return []


def initialize_header(hdr=None):
    """
    Initialize a FITS header.

    Args:
        hdr (`astropy.io.fits.Header`, optional):
            Header object to update with basic summary
            information. The object is modified in-place and also
            returned. If None, an empty header is instantiated,
            edited, and returned.

    Returns:
        `astropy.io.fits.Header`: The initialized (or edited)
        fits header.
    """
    # Add versioning; this hits the highlights but should it add
    # the versions of all packages included in the requirements.txt
    # file?
    if hdr is None:
        hdr = fits.Header()
    hdr['VERSPYT'] = ('.'.join([ str(v) for v in sys.version_info[:3]]), 'Python version')
    hdr['VERSNPY'] = (numpy.__version__, 'Numpy version')
    hdr['VERSSCI'] = (scipy.__version__, 'Scipy version')
    hdr['VERSAST'] = (astropy.__version__, 'Astropy version')
    hdr['VERSSKL'] = (sklearn.__version__, 'Scikit-learn version')
    hdr['VERSPYP'] = (pypeit.__version__, 'PypeIt version')

    # Save the date of the reduction
    hdr['DATE'] = (time.strftime('%Y-%m-%d',time.gmtime()), 'UTC date created')

    # TODO: Anything else?

    # Return
    return hdr

def header_version_check(hdr, warning_only=True):
    """
    Check the package versions in the header match the system versions.

    .. note::
        The header must contain the keywords written by
        :func:`initialize_header`.

    Args:
        hdr (`astropy.io.fits.Header`):
            The header to check
        warning_only (:obj:`bool`, optional):
            If the versions are discrepant, only throw a warning
            instead of raising an exception.

    Returns:
        :obj:`bool`: Returns True if the check was successful, False
        otherwise. If `warning_only` is False, the method will either
        raise an exception or return True.

    Raises:
        ValueError:
            Raised if `warning_only` is False and the system versions
            are different from those logged in the header.
    """
    # Compile the packages and versions to check
    packages = ['python', 'numpy', 'scipy', 'astropy', 'sklearn', 'pypeit']
    hdr_versions = [hdr['VERSPYT'], hdr['VERSNPY'], hdr['VERSSCI'], hdr['VERSAST'], hdr['VERSSKL'],
                    hdr['VERSPYP']]
    sys_versions = ['.'.join([ str(v) for v in sys.version_info[:3]]), numpy.__version__,
                    scipy.__version__, astropy.__version__, sklearn.__version__,
                    pypeit.__version__]

    # Run the check and either issue warnings or exceptions
    all_identical = True
    for package, hdr_version, sys_version in zip(packages, hdr_versions, sys_versions):
        if version.parse(hdr_version) != version.parse(sys_version):
            all_identical = False
            msg = '{0} version used to create the file ({1}) '.format(package, hdr_version) \
                        + 'does not match the current system version ({0})!'.format(sys_version)
            if warning_only:
                # TODO: I had to change pypeit/__init__.py to get these
                # to show up. We eventually need to make pypmsgs play
                # nice with warnings and other logging, or just give up
                # on pypmsgs...
                warnings.warn(msg)
            else:
                raise ValueError(msg)
    # Return if all versions are identical
    return all_identical

