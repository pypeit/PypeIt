"""
Implements classes and function for the PypeIt data model.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import warnings

from IPython import embed

import numpy as np

from astropy.io import fits
from astropy.table import Table

from pypeit import io

class DataContainer:
    """
    Defines an abstract class for holding and manipulating data.

    The primary utilities of the class are:
        - Attributes can be accessed normally or as expected for a :obj:`dict`
        - Attributes and items are restricted to conform to a specified data model.

    This abstract class should only be used as a base class.

    Derived classes must do the following:

        - Define a datamodel
        - Provide an :func:`__init__` method that defines the
          instantiation calling sequence and passes the relavant
          dictionary to this base-class instantiation.
        - Provide a :func:`_validate` method, if necessary, that
          processes the data provided in the `__init__` into a
          complete instantiation of the object. This method and the
          :func:`__init__` method are the *only* places where
          attributes can be added to the class.
        - Provide a :func:`_bundle` method that reorganizes the
          datamodel into partitions that can be written to one or
          more fits extensions. More details are provided in the
          description of :func:`_bundle`.
        - Provide a :func:`_parse` method that parses information in
          one or more fits extensions into the appropriate datamodel.
          More details are provided in the description of
          :func:`_parse`.

    Args:
        d (:obj:`dict`):
            Dictionary to copy to the internal attribute dictionary.
    """

    # Define the data model
    datamodel = None
    """
    Provides the class data model. This is undefined in the abstract
    class and should be overwritten in the derived classes.

    .. note::
        Data model elements should *not* be dicts themselves.
    """
    def __init__(self, d):
        # Data model must be defined
        if self.datamodel is None:
            raise ValueError('Data model for {0} is undefined!'.format(self.__class__.__name__))

        # Input dictionary cannot have keys that do not exist in the
        # data model
        if not np.all(np.isin(list(d.keys()), list(self.datamodel.keys()), assume_unique=True)):
            raise AttributeError('Coding error: Initialization arguments do not match data model!')

        # Copy the dictionary
#        self.__dict__ = copy.deepcopy(d)
        # Ensure the dictionary has all the expected keys
        self.__dict__.update(dict.fromkeys(self.datamodel.keys()))
        # Assign the values provided by the input dictionary
        self.__dict__.update(d)

        # Validate the object
        # TODO: _validate isn't the greatest name for this method...
        self._validate()

        # TODO: Confirm elements have otype and atypes consistent with
        # the data model?

        # Finalize the instantiation.
        # NOTE: The key added to `__dict__` by this call is always
        # `_DataContainer__initialised`, regardless of whether or not
        # the call to this `__init__` is from the derived class.
        # Originally, I had been checking if the result from the
        # `_init_key` method below was in `__dict__`, but this gives
        # different names for the derived classes. In the end, I check
        # if `_DataContainer__initialised` is in `self.__dict__` (as
        # done in the SpecObj example). Is there a way we could just
        # check a boolean instead?
        self.__initialised = True

    def _validate(self):
        """
        Validate the data container.

        The purpose of this function is to check the input data
        provided by the instantiation and flesh out any details of
        the object.

        Derived classes should override this function, unless there
        is nothing to validate.

        Attributes can be added to the object in this function
        because it is called before the datamodel is frozen.
        """
        pass

    def _bundle(self, ext=None, transpose_arrays=False):
        """
        Bundle the data into a series of objects to be written to
        fits HDU extensions.

        The returned object must be a list. The list items should be
        one of the following:

            - a dictionary with a single key/item pair, where the key
              is the name for the extension and the item is the data
              to be written.

            - a single item (object) to be written, which will be
              written in the provided order without any extension
              names (although see the caveat in
              :func:`pypeit.io.dict_to_hdu` for dictionaries with
              single array or `astropy.table.Table`_ items).

        The item to be written can be a single array for an ImageHDU,
        an `astropy.table.Table`_ for a BinTableHDU, or a dictionary;
        see :func:`pypeit.io.write_to_hdu`.

        For how these objects parsed into the HDUs, see
        :func:`to_hdu`.

        The default behavior implemented by this base class just
        parses the attributes into a single dictionary based on the
        datamodel and is returned such that it all is written to a
        single fits extension.  Note that this will fault if the datamodel contains:

            - a dictionary object
            - more than one `astropy.table.Table`_,
            - an `astropy.table.Table`_ and an array-like object
              (:obj:`list` or `numpy.ndarray`_), or

        Certain **restrictions** apply to how the data can be bundled
        for the general parser implementation (:func:`_parse`) to
        work correctly. These restriction are:

            - The shape and orientation of any input arrays are
              assumed to be correct.

            - Datamodel keys for arrays or `astropy.table.Table`_
              objects written to an HDU should match the HDU
              extension name. Otherwise, the set of HDU extension
              names and datamodel keys **must** be unique.

            - Datamodel keys will be matched to header values:
              :func:`to_hdu` will write any items in a dictionary
              that is an integer, float, or string (specific numpy
              types or otherwise) to the header. This means header
              keys in **all** extensions should be unique and should
              not be the same as any extension name.

            - Datamodels can contain tuples, but these must be
              reasonably converted to strings such they are included
              in (one of) the HDU header(s).

        Args:
            ext (:obj:`str`, optional):
                Name for the HDU extension. If None, no name is
                given.
            transpose_arrays (:obj:`bool`, optional):
                Transpose the arrays before writing them to the HDU.
                This option is mostly meant to correct orientation of
                arrays meant for tables so that the number of rows
                match.
            
        Returns:
            :obj:`list`: A list of dictionaries, each list element is
            written to its own fits extension. See the description
            above.

        Raises:
            TypeError:
                Raised if the provided ``ext`` is not a string.
        """
        d = {}
        for key in self.keys():
            if transpose_arrays and self.datamodel[key]['otype'] == np.ndarray:
                d[key] = self[key].T
            elif self.datamodel[key]['otype'] == tuple:
                d[key] = str(self[key])
            else:
                d[key] = self[key]
        return [d] if ext is None else [{ext:d}]

    @classmethod
    def _parse(cls, hdu, ext=None, transpose_table_arrays=False):
        """
        Parse data read from a set of HDUs.

        This method is the counter-part to :func:`_bundle`, and
        parses data from the HDUs into a dictionary that can be used
        to instantiate the object.

        .. warning::

            - Beware that this may read data from the provided HDUs
              that could then be changed by :func:`_validate`.
              Construct your :func:`_validate` methods carefully!

            - Although :func:`_bundle` methods will likely need to be
              written for each derived class, this parsing method is
              very general to what :class:`DataContainer` can do.
              Before overwriting this function in a derived class,
              make sure and/or test that this method doesn't meet
              your needs, and then tread carefully regardless.

            - Because the `astropy.table.Table`_ methods are used
              directly, any metadata associated with the Table will
              also be included in the HDUs construced by
              :func:`to_hdu`. However, the
              `astropy.table.Table.read`_ method always returns the
              metadata with capitalized keys. This means that,
              regardless of the capitalization of the metadata
              keywords when the data is written, **they will be
              upper-case when read by this function**!

        .. note::

            Currently, the only test for the object type (``otype``)
            given explicitly by the class datamodel is to check if
            the type should be a tuple. If it is, the parser reads
            the string from the header and evaluates it so that it's
            converted to a tuple on output. See the restrictions
            listed for :func:`_bundle`.

            All the other type conversions are implicit or based on
            the HDU type.

        Args:
            hdu (`astropy.io.fits.HDUList`_, `astropy.io.fits.ImageHDU`_,
                 `astropy.io.fits.BinTableHDU`_):
                The HDU(s) to parse into the instantiation dictionary.
            ext (:obj:`int`, :obj:`str`, :obj:`list`, optional):
                One or more extensions with the data. If None, the
                function trolls through the HDU(s) and parses the
                data into the datamodel.
            transpose_table_arrays (:obj:`bool`, optional):
                Tranpose *all* the arrays read from any binary
                tables. This is meant to invert the use of
                ``transpose_arrays`` in :func:`_bound`.

        Returns:
            :obj:`dict`: Dictionary used to instantiate the object.

        Raises:
            TypeError:
                Raised if ``ext``, or any of its elements if it's a
                :obj:`list`, are not either strings or integers.
        """
        # Get the (list of) extension(s) to parse
        if ext is not None:
            if not isinstance(ext, (str, int, list)):
                raise TypeError('Provided ext object must be a str, int, or list.')
            if isinstance(ext, list):
                for e in ext:
                    if not isinstance(ext, (str, int)):
                        raise TypeError('Provided ext elements  must be a str or int.')
        if ext is None and isinstance(hdu, fits.HDUList):
            ext = [h.name for h in hdu]

        # Allow user to provide single HDU
        if isinstance(hdu, (fits.ImageHDU, fits.BinTableHDU)):
            if ext is not None:
                warnings.warn('Only one HDU provided; extension number/name is irrelevant.')
            ext = [0]
            _hdu = [hdu]
        else:
            _hdu = hdu

        _ext = np.atleast_1d(ext)

        # Construct instantiation dictionary
        d = dict.fromkeys(cls.datamodel.keys())

        # NOTE: The extension and keyword comparisons are complicated
        # because the fits standard is to force these all to be
        # capitalized, while the datamodel doesn't (currently)
        # implement this restriction.

        # HDUs can have dictionary elements directly.
        keys = np.array(list(d.keys()))
        indx = np.isin([key.upper() for key in keys], _ext)
        if np.any(indx):
            for e in keys[indx]:
                d[e] = hdu[e.upper()].data if isinstance(hdu[e.upper()], fits.ImageHDU) \
                        else Table.read(hdu[e.upper()])

        for e in _ext:
            # Check for header elements
            indx = np.isin([key.upper() for key in keys], list(hdu[e].header.keys()))
            if np.any(indx):
                for key in keys[indx]:
                    d[key] = hdu[e].header[key.upper()] if cls.datamodel[key]['otype'] != tuple \
                                else eval(hdu[e].header[key.upper()])
            # Parse BinTableHDUs
            if isinstance(hdu[e], fits.BinTableHDU):
                # If the length of the table is 1, assume the table
                # data had to be saved as a single row because of shape
                # differences.
                single_row = len(hdu[e].data) == 1
                for key in hdu[e].columns.names:
                    if key in cls.datamodel.keys():
                        d[key] = hdu[e].data[key][0] if single_row else hdu[e].data[key]
                        if transpose_table_arrays:
                            d[key] = d[key].T
        return d
                
#    @classmethod
#    def _init_key(cls):
#        """
#        The key in the internal dict that establishes that the object
#        is initialized.
#        """
#        return '_{0}__initialised'.format(cls.__name__)

    def __getattr__(self, item):
        """Maps values to attributes.
        Only called if there *isn't* an attribute with this name
        """
        try:
            return self.__getitem__(item)
        except KeyError:
            raise AttributeError('{0} is not an attribute of {1}!'.format(item,
                                    self.__class__.__name__))

    def __setattr__(self, item, value):
        """
        Set the attribute value.

        Attributes are restricted to those defined by the datamodel.
        """
        # TODO: It seems like it would be faster to check a boolean
        # attribute. Is that possible?
        if '_DataContainer__initialised' not in self.__dict__:
            # Allow new attributes to be added before object is
            # initialized
            dict.__setattr__(self, item, value)
        else:
            # Otherwise, set as an item
            try:
                self.__setitem__(item, value)
            except KeyError as e:
                # Raise attribute error instead of key error
                raise AttributeError(e)

#        # TODO: It seems like it would be faster to check a boolean
#        # attribute. Is that possible?
#        if '_DataContainer__initialised' not in self.__dict__:
#            # Allow new attributes to be added before object is
#            # initialized
#            dict.__setattr__(self, item, value)
#        elif item in self.__dict__:
#            # Only set attributes that already exist
#            dict.__setattr__(self, item, value)
#        else:
#            # Otherwise, set as an item
#            self.__setitem__(item, value)

    def __setitem__(self, item, value):
        """
        Access and set an attribute identically to a dictionary item.

        Items are restricted to those defined by the datamodel.
        """
        if item not in self.datamodel.keys():
            raise KeyError('Key {0} not part of the data model'.format(item))
        if not isinstance(value, self.datamodel[item]['otype']):
            raise TypeError('Incorrect data type for {0}!'.format(item) + 
                            'Allowed type(s) are: {0}'.format(self.datamodel[item]['otype']))
        self.__dict__[item] = value

    def __getitem__(self, item):
        return self.__dict__[item]

    def keys(self):
        """
        Return the keys for the data objects.

        Returns:
            :obj:`dict_keys`: The iterable with the data model keys.
        """
        return self.datamodel.keys()

    def to_hdu(self, add_primary=False):
        """
        Construct one or more HDU extensions with the data.

        The organization of the data into separate HDUs is performed
        by :func:`_bundle`, which returns a list of objects. The type
        of each element in the list determines how it is handled. If
        the object is a dictionary with a single key/item pair, the
        key is used as the extension header. Otherwise, the objects
        are written to unnamed extensions. The object or dictionary
        item is passed to :func:`pypeit.io.write_to_hdu` to construct
        the HDU.

        Args:
            add_primary (:obj:`bool`, optional):
                If False, the returned object is a simple
                :obj:`list`, with a list of HDU objects (either
                `astropy.io.fits.ImageHDU`_ or
                `astropy.io.fits.BinTableHDU`_). If true, the method
                constructs an `astropy.io.fits.HDUList` with a
                primary HDU, such that this call::

                    hdu = fits.HDUList([fits.PrimaryHDU()] + self.to_hdu())

                and this call::

                    hdu = self.to_hdu(add_primary=True)

                are identical.

        Returns:
            :obj:`list`, `astropy.io.fits.HDUList`_: A list of HDUs,
            where the type depends on the value of ``add_primary``.
        """
        # Bundle the data
        data = self._bundle()

        # Construct the list of HDUs
        hdu = []
        for d in data:
            if isinstance(d, dict) and len(d) == 1:
                ext = list(d.keys())[0]
                hdu += [io.write_to_hdu(d[ext], name=ext)]
            else:
                hdu += [io.write_to_hdu(d)]
        return fits.HDUList([fits.PrimaryHDU()] + hdu) if add_primary else hdu

    @classmethod
    def from_hdu(cls, hdu):
        """
        Instantiate the object from an HDU extension.

        This is primarily a wrapper for :func:`_parse`.

        Args:
            hdu (`astropy.io.fits.HDUList`_, `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_):
                The HDU(s) with the data to use for instantiation.
        """
        # NOTE: We can't use `cls(cls._parse(hdu))` here because this
        # will call the `__init__` method of the derived class and we
        # need to use the `__init__` of the base class instead. So
        # below, I get an empty instance of the derived class using
        # `__new__`, call the parent `__init__`, and then return the
        # result. This is the first time I've used `__new__` so we may
        # want to tread carefully, but it seems to work.
        self = super().__new__(cls)
        super(cls, self).__init__(cls._parse(hdu))
        return self

    def to_file(self, ofile, overwrite=False, checksum=True):
        """
        Write the data to a file.

        This is a convenience wrapper for :func:`to_hdu` and
        :func:`pypeit.io.write_to_fits`. The output is always placed
        in the 2nd extension; the first (primary) extension is always
        empty.

        Args:
            ofile (:obj:`str`):
                Fits file for the data. File names with '.gz'
                extensions will be gzipped; see
                :func:`pypeit.io.write_to_fits`.
            overwrite (:obj:`bool`, optional):
                Flag to overwrite any existing file.
            checksum (:obj:`bool`, optional):
                Passed to `astropy.io.fits.HDUList.writeto`_ to add
                the DATASUM and CHECKSUM keywords fits header(s).
        """
        io.write_to_fits(self.to_hdu(add_primary=True), ofile, overwrite=overwrite,
                         checksum=checksum)

    @classmethod
    def from_file(cls, ifile):
        """
        Instantiate the object from an extension in the specified fits file.

        This is a convenience wrapper for :func:`from_hdu`.
        
        Args:
            ifile (:obj:`str`):
                Fits file with the data to read

        Raises:
            FileNotFoundError:
                Raised if the specified file does not exist.
        """
        if not os.path.isfile(ifile):
            raise FileNotFoundError('{0} does not exist!'.format(ifile))
        with fits.open(ifile) as hdu:
            return cls.from_hdu(hdu)

