# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides a set of I/O routines.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import os
import sys
import warnings
import gzip
import shutil
from packaging import version

import numpy

from configobj import ConfigObj

from astropy.io import fits
from astropy.table import Table

from pypeit import par, msgs

# These imports are largely just to make the versions available for
# writing to the header. See `initialize_header`
import scipy
import astropy
import sklearn
import pypeit
import time


from IPython import embed

# TODO -- Move this module to core/

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

# TODO: Should probably rename this since it's no longer only used for
# record arrays.
def rec_to_fits_type(col_element, single_row=False):
    """
    Return the string representation of a fits binary table data type
    based on the provided column element.

    Args:
        col_element (`numpy.ndarray`_):
            The example data to write to a
            `astropy.io.fits.BinTableHDU`_ used to determine the
            column format.
        single_row (:obj:`bool`, optional):
            Flag that the provided object is the data written to a
            single row for the `astropy.io.fits.BinTableHDU`_ column.

    Returns:
        str: String representation of the format for the column.
    """
    _col_element = col_element if single_row else col_element[0]
    n = 1 if len(_col_element.shape) == 0 else _col_element.size
    if col_element.dtype.type in [bool, numpy.bool, numpy.bool_]:
        return '{0}L'.format(n)
    if col_element.dtype.type == numpy.uint8:
        return '{0}B'.format(n)
    if col_element.dtype.type in [numpy.int16, numpy.uint16]:
        return '{0}I'.format(n)
    if col_element.dtype.type in [numpy.int32, numpy.uint32]:
        return '{0}J'.format(n)
    if col_element.dtype.type in [numpy.int64, numpy.uint64]:
        return '{0}K'.format(n)
    if col_element.dtype.type == numpy.float32:
        return '{0}E'.format(n)
    if col_element.dtype.type == numpy.float64:
        return '{0}D'.format(n)
    if col_element.dtype.name == 'float32':  # JXP -- Hack for when slit edges are modified in the Flat making
        return '{0}E'.format(n)
    if col_element.dtype.name == 'float64':  # JXP -- Hack for when slit edges are modified in the Flat making
        return '{0}D'.format(n)

    # If it makes it here, assume its a string
    l = int(col_element.dtype.str[col_element.dtype.str.find('U')+1:])
#    return '{0}A'.format(l) if n==1 else '{0}A{1}'.format(l*n,l)
    return '{0}A'.format(l*n)


def rec_to_fits_col_dim(col_element, single_row=False):
    """
    Return the string representation of the dimensions for the fits
    table column based on the provided column element.

    The shape is inverted because the first element is supposed to be
    the most rapidly varying; i.e. the shape is supposed to be written
    as row-major, as opposed to the native column-major order in python.

    Args:
        col_element (`numpy.ndarray`_):
            The example data to write to a
            `astropy.io.fits.BinTableHDU`_ used to determine the
            column dimension.
        single_row (:obj:`bool`, optional):
            Flag that the provided object is the data written to a
            single row for the `astropy.io.fits.BinTableHDU`_ column.

    Returns:
        str: String representation of the column dimensions. Return
        None if the object is not multidimensional.
    """
    _col_element = col_element if single_row else col_element[0]
    return None if len(_col_element.shape) < 2 else str(_col_element.shape[::-1])


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


def initialize_header(hdr=None, primary=False):
    """
    Initialize a FITS header.

    Args:
        hdr (`astropy.io.fits.Header`, optional):
            Header object to update with basic summary
            information. The object is modified in-place and also
            returned. If None, an empty header is instantiated,
            edited, and returned.
        primary (bool, optional):
            If true and hdr is None, generate a Primary header

    Returns:
        `astropy.io.fits.Header`: The initialized (or edited)
        fits header.
    """
    # Add versioning; this hits the highlights but should it add
    # the versions of all packages included in the requirements.txt
    # file?
    if hdr is None:
        if primary:
            hdr = fits.PrimaryHDU().header
        else:
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


def dict_to_hdu(d, name=None, hdr=None, force_to_bintbl=False):
    """
    Write a dictionary to a fits HDU.

    Elements in the dictionary that are integers, floats, or strings
    (specific numpy types or otherwise) are written to the HDU
    header. The header keywords are identical to the dictionary keys.

    If any of the elements in the dictionary are an
    `astropy.table.Table`_, that dictionary can *only* contain that
    table and single values that will be written to the extension
    header. That is, there can be only one `astropy.table.Table`_
    element, and none of the elements can be a :obj:`list` or
    `numpy.ndarray`_. By default the extension name is the dictionary
    key for the `astropy.table.Table`_ item; this can be overridden
    using the ``name`` argument.

    Elements in the dictionary that are a list or a `numpy.ndarray`_
    are written as either an image (if there is only one array and a
    binary table is not specifically requested using
    ``force_to_bintbl``) or a series of table columns. The lists are
    assumed to be interpretable as the ``array`` argument of
    `astropy.io.fits.Column`_ (for a table) or the ``data`` argument
    of `astropy.io.fits.ImageHDU`_ (for an image).

        - If an image is to be written, the extension name, by
          default, is the dictionary key for the array item; this can
          be overridden using the ``name`` argument.

        - If a table is to be written, the method checks that the
          relevant arrays have a consistent number of rows. If they
          do not, the format and dimensions of the table written are
          set so that the arrays are contained in a single row. The
          column names in the table are identical to the dictionary
          keywords. In this case, ``name`` must be provided if you
          want the extension to have a name; there is no default
          name.

    Args:
        d (:obj:`dict`):
            Dictionary with data to write to the
            `astropy.io.fits.BinTableHDU`_.
        name (:obj:`str`, optional):
            Name to give the HDU extension. If None and the input is
            a dictionary with a single array or
            `astropy.table.Table`_ to write, the name of the
            extension is the relevant dictionary keyword. Any
            provided value for ``name`` will override this behavior.
            If the provided dictionary is used to construct a table,
            where the dictionary keys are used for the table column
            names, there is no default name for the extension (i.e.,
            no extension name is used if ``name is None``).
        hdr (`astropy.io.fits.Header`_, optional):
            Base-level header to include in the HDU. If None, an
            empty header is used and then added to.
        force_to_bintbl (:obj:`bool`, optional):
            Force a BinTableHDU to be constructed instead of an
            ImageHDU when either there are no arrays or tables to
            write or only a single array is provided.

    Returns:
        `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_:
        HDU with the data. An `astropy.io.fits.ImageHDU`_ object is
        returned if there is 1 (or fewer) array-like objects in the
        dictionary. Otherwise, an `astropy.io.fits.BinTableHDU`_
        object is returned with the data.

    Raises:
        TypeError:
            Raised if the input object is not a dictionary or the
            method cannot interpret how to use an element of the
            dictionary.
        ValueError:
            Raised if dictionary contains another dictionary, more
            than one `astropy.table.Table`_ object, or both an
            `astropy.table.Table`_ and an array-like object
            (:obj:`list` or `numpy.ndarray`_).
    """
    # Check the input is a dictionary (not very pythonic...)
    if not isinstance(d, dict):
        raise TypeError('Input must be a dictionary.')
    # Check the dictionary contents
    ndict = numpy.sum([isinstance(d[key], dict) for key in d.keys()])
    if ndict > 0:
        raise ValueError('Cannot write nested dictionaries.')
    ntab = numpy.sum([isinstance(d[key], Table) for key in d.keys()])
    if ntab > 1:
        raise ValueError('Cannot write dictionaries with more than one astropy.table.Table.')
    narr = numpy.sum([isinstance(d[key], (list, numpy.ndarray)) for key in d.keys()])
    if ntab > 0 and narr > 0:
        raise ValueError('Cannot write dictionaries with both arrays and Tables.')

    # Write any header data and find arrays and Tables
    _hdr = fits.Header() if hdr is None else hdr.copy()
    array_keys = []
    table_keys = []
    for key in d.keys():
        if d[key] is None:
            continue
        # TODO: may be better to do this
        #   isinstance(d[key], (collections.Sequence, numpy.ndarray)):
        # This ignores the defined otype...
        if isinstance(d[key], (list, numpy.ndarray)):
            array_keys += [key]
        elif isinstance(d[key], Table):
            table_keys += [key]
        elif isinstance(d[key], (int, numpy.integer, float, numpy.floating, str)):
            _hdr[key.upper()] = d[key]
        else:
            raise TypeError('Do not know how to write object with type {0}'.format(type(d[key])))

    # If there aren't any arrays or tables, return an empty ImageHDU or
    # BinTableHDU with just the header data.
    if len(array_keys) < 1 and len(table_keys) < 1:
        return fits.BinTableHDU(header=_hdr, name=name) if force_to_bintbl \
                    else fits.ImageHDU(header=_hdr, name=name)

    # If there's only a single array, return it in an ImageHDU or, if
    # requested, a BinTableHDU
    if len(array_keys) == 1 and not force_to_bintbl:
        return fits.ImageHDU(data=d[array_keys[0]], header=_hdr,
                             name=array_keys[0] if name is None else name)

    # If there's only a single Table, return it in a BinTableHDU
    if len(table_keys) == 1:
        # TODO: If we pass hdr directly, does this call include any
        # table meta?
        return fits.BinTableHDU(data=d[table_keys[0]], header=_hdr,
                                name=table_keys[0] if name is None else name)

    # Only remaining option is to build a BinTableHDU based on the
    # dictionary contents.

    # Do all arrays have the same number of rows?
    single_row = len(numpy.unique([numpy.asarray(d[key]).shape[0] for key in array_keys])) > 1

    # If the number of rows is inconsistent, save the data in a single
    # row. Otherwise, save the data as a multi-row table.
    cols = []
    for key in array_keys:
        _d = numpy.asarray(d[key])
        cols += [fits.Column(name=key,
                             format=rec_to_fits_type(_d, single_row=single_row),
                             dim=rec_to_fits_col_dim(_d, single_row=single_row),
                             array=numpy.expand_dims(_d, 0) if single_row else _d)]
    return fits.BinTableHDU.from_columns(cols, header=_hdr, name=name)


def read_spec2d_file(ifile, filetype='coadd2d'):
    """
    Read a PypeIt file of type "filetype", akin to a standard PypeIt file.

    .. todo::

        - Need a better description of this.  Probably for the PypeIt
          file itself, too!

    The top is a config block that sets ParSet parameters.  The
    spectrograph is required.

    Args:
        ifile (:obj:`str`):
          Name of the config file
        filetype (:obj:`str`):
          Type of config file being read (e.g. coadd2d, coadd3d).
          This is only used for printing messages.

    Returns:
        tuple: Returns three objects: (1) The name of the spectrograph as
        a string, (2) the list of configuration lines used to modify the
        default :class`pypeit.par.pypeitpar.PypeItPar` parameters, and
        (3) the list of spec2d files to combine.
    """

    # Read in the pypeit reduction file
    msgs.info('Loading the {0:s} file'.format(filetype))
    lines = par.util._read_pypeit_file_lines(ifile)
    is_config = numpy.ones(len(lines), dtype=bool)

    # Parse the spec2d block
    spec2d_files = []
    s, e = par.util._find_pypeit_block(lines, 'spec2d')
    if s >= 0 and e < 0:
        msgs.error("Missing 'spec2d end' in {0}".format(ifile))
    for line in lines[s:e]:
        prs = line.split(' ')
        # TODO: This needs to allow for the science directory to be
        # defined by the user.
        #spec2d_files.append(os.path.join(os.path.basename(prs[0])))
        spec2d_files.append(prs[0])
    is_config[s-1:e+1] = False
    # Construct config to get spectrograph
    cfg_lines = list(lines[is_config])
    cfg = ConfigObj(cfg_lines)
    spectrograph_name = cfg['rdx']['spectrograph']

    # Return
    return spectrograph_name, cfg_lines, spec2d_files


def write_to_hdu(d, name=None, hdr=None, force_to_bintbl=False):
    """
    Write the input to an astropy.io.fits HDU extension.

    List and numpy.ndarray items are written as an ImageHDU,
    `astropy.table.Table`_ items are written as a BinTableHDU, and
    dictionaries are passed to :func:`dict_to_hdu`.
    
    .. warning::

        If ``d`` is a list, the method assumes that the list is
        essentially an array that can be sensibly converted to a
        `numpy.ndarray`_.

    Args:
        d (:obj:`dict`, :obj:`list`, `numpy.ndarray`_, `astropy.table.Table`_):
            Object to write to the HDU.
        name (:obj:`str`, optional):
            Name for the HDU extension.
        hdr (`astropy.io.fits.Header`_, optional):
            Header to include in the HDU.
        force_to_bintbl (bool, optional):
            Force dict into a BinTableHDU

    Returns:
        `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_:
        HDU with the data.

    Raises:
        TypeError:
            Raised if the input object is not one of the allowed
            types.

    """
    if isinstance(d, dict):
        return dict_to_hdu(d, name=name, hdr=hdr, force_to_bintbl=force_to_bintbl)
    if isinstance(d, Table):
        return fits.BinTableHDU(data=d, name=name, header=hdr)
    if isinstance(d, (numpy.ndarray, list)):
        return fits.ImageHDU(data=numpy.asarray(d), name=name, header=hdr)
    raise TypeError('Input must be a dictionary, astropy.table.Table, list, or numpy.ndarray.')


def write_to_fits(d, ofile, name=None, hdr=None, overwrite=False, checksum=True):
    """
    Write the provided object to a fits file.

    This is either a convenience wrapper for :func:`write_to_hdu`
    that adds a primary HDU and writes the result to the provided
    file, or a convenience wrapper for an already formed
    `astropy.io.fits.HDUList`_ passed as (``d``).

    If the provided file name includes the '.gz' extension, the file
    is first written using `astropy.io.fits.HDUList.writeto`_ and
    then compressed using :func:`compress_file`.
    
    .. note::

        - If the root directory of the output does *not* exist, this
          method will create it.
        - Compressing the file is generally slow, but following the
          two-step process of running
          `astropy.io.fits.HDUList.writeto`_ and then
          :func:`compress_file` is generally faster than having
          `astropy.io.fits.HDUList.writeto`_ do the compression,
          particularly for files with many extensions (or at least
          this was true in the past).

    Args:
        d (:obj:`dict`, :obj:`list`, `numpy.ndarray`_, `astropy.table.Table`_, `astropy.io.fits.HDUList`_):
            Object to write to the HDU. See :func:`write_to_hdu`.
        ofile (:obj:`str`):
            File name (path) for the fits file.
        name (:obj:`str`, optional):
            Name for the extension with the data. If None, the
            extension is not given a name. However, if the input
            object is a dictionary, see :func:`dict_to_hdu` for how
            the name will overwrite any dictionary keyword associated
            with the data to write.
        hdr (`astropy.io.fits.Header`_, optional):
            Base-level header to use for *all* HDUs.
        overwrite (:obj:`bool`, optional):
            Overwrite any existing file.
        checksum (:obj:`bool`, optional):
            Passed to `astropy.io.fits.HDUList.writeto`_ to add the
            DATASUM and CHECKSUM keywords fits header(s).
    """
    if os.path.isfile(ofile) and not overwrite:
        raise FileExistsError('File already exists; to overwrite, set overwrite=True.')
    
    root = os.path.split(os.path.abspath(ofile))[0]
    if not os.path.isdir(root):
        warnings.warn('Making root directory for output file: {0}'.format(root))
        os.makedirs(root)

    # Determine if the file should be compressed
    _ofile = ofile[:ofile.rfind('.')] if ofile.split('.')[-1] == 'gz' else ofile

    _hdr = initialize_header() if hdr is None else hdr.copy()

    # Construct the hdus and write the fits file.
    fits.HDUList(d if isinstance(d, fits.HDUList) else 
                 [fits.PrimaryHDU(header=_hdr)] + [write_to_hdu(d, name=name, hdr=_hdr)]
                 ).writeto(_ofile, overwrite=True, checksum=checksum)

    # Compress the file if the output filename has a '.gz' extension;
    # this is slow but still faster than if you have astropy.io.fits do
    # it directly
    # TODO: use pypmsgs?
    if _ofile is not ofile:
        pypeit.msgs.info('Compressing file: {0}'.format(_ofile))
        compress_file(_ofile, overwrite=True)
    pypeit.msgs.info('File written to: {0}'.format(ofile))


def hdu_iter_by_ext(hdu, ext=None, hdu_prefix=None):
    """
    Convert the input to lists that can be iterated through by an
    extension index/name.

    If ``hdu`` is an `astropy.io.fits.HDUList`_ on input, it is
    simply returned; otherwise, the 2nd returned item is a
    single-element :obj:`list` with the provided HDU.

    If ``ext`` is None and ``hdu`` is not an
    `astropy.io.fits.HDUList`_, the returned list just selects the
    individual HDU provided (i.e., ``ext = [0]``). If ``ext`` is None
    and ``hdu`` *is* an `astropy.io.fits.HDUList`_, the returned list
    of extensions includes all extensions in the provided ``hdu``.

    .. warning::

        The method does not check that all input ``ext`` are valid
        for the provided ``hdu``!

    Args:
        hdu (`astropy.io.fits.HDUList`_, `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_):
            The HDU(s) to iterate through.
        ext (:obj:`int`, :obj:`str`, :obj:`list`, optional):
            One or more extensions to include in the iteration. If
            None, the returned list will enable iteration through all
            HDU extensions.
        hdu_prefix (:obj:`str`, optional):
            In addition to the restricted list of extensions
            (``ext``), only include extensions with this prefix.

    Returns:
        :obj:`tuple`: Returns two objects: a :obj:`list` with the
        extensions to iterate through and either a :obj:`list` or an
        `astropy.io.fits.HDUList`_ with the list of HDUs.

    Raises:
        TypeError:
            Raised if ``ext`` is not a string, integer, or list, if
            any element of ``ext`` is not a string or integer, or if
            ``hdu`` is not one of the approved types.
    """
    if not isinstance(hdu, (fits.HDUList, fits.ImageHDU, fits.BinTableHDU)):
        raise TypeError('Provided hdu has incorrect type: {0}'.format(type(hdu)))
    # Check that the input ext has valid types
    if ext is not None:
        if not isinstance(ext, (str, int, list)):
            raise TypeError('Provided ext object must be a str, int, or list.')
        if isinstance(ext, list):
            for e in ext:
                if not isinstance(ext, (str, int)):
                    raise TypeError('Provided ext elements  must be a str or int.')
    if ext is None and isinstance(hdu, fits.HDUList):
        ext = [h.name if h.name != '' else i for i,h in enumerate(hdu)]

    # Further restrict ext to only include those with the designated
    # prefix
    if hdu_prefix is not None:
        if isinstance(hdu, fits.HDUList):
            ext = [e for e in ext if not isinstance(e, (int, numpy.integer))
                                 and e.startswith(hdu_prefix)]

    # Allow user to provide single HDU
    if isinstance(hdu, (fits.ImageHDU, fits.BinTableHDU)):
        ext = [0]
        _hdu = [hdu]
        if hdu_prefix is not None:
            if hdu_prefix not in hdu.name:
                raise ValueError("Bad hdu_prefix for this HDU!")
    else:
        _hdu = hdu

    return ext if isinstance(ext, list) else [ext], _hdu

