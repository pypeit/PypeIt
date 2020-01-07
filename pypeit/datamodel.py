"""
Implements classes and function for the PypeIt data model.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import copy

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit import io

class DataContainer:
    """
    Defines an abstract class for holding and manipulating data.

    The primary utilities of the class are:
        - Attributes can be accessed normally or as expected for a :obj:`dict`
        - Attributes and items are restricted to conform to a specified data model.

    This abstract class should only be used as a base class.

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
        self.__dict__ = copy.deepcopy(d)

        # Validate the object
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

        Derived classes should override this function, unless there
        is nothing to validate.
        """
        pass

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
        # TODO: Not sure this is needed in this implementation
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
        elif item in self.__dict__:
            # Only set attributes that already exist
            dict.__setattr__(self, item, value)
        else:
            # Otherwise, set as an item
            self.__setitem__(item, value)

    def __setitem__(self, item, value):
        """
        Access and set an attribute identically to a dictionary item.

        Items are restricted to those defined by the datamodel.
        """
        if item not in self.datamodel.keys():
            raise AttributeError('Attribute {0} not part of the data model'.format(item))
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

    def to_hdu(self, name=None):
        """
        Construct an `astropy.io.fits.BinTableHDU`_ with the data.

        The function uses the datamodel to construct an
        `astropy.io.fits.BinTableHDU`_ with a single row that
        contains the data:

            - Any array-like objects are written as columns in the
              table, with the column names set by the datamodel keys.

            - Objects with individual values (int, float, str, etc.)
            are written to the extension header.

        Args:
            name (:obj:`str`, optional):
                Extension name given to the
                `astropy.io.fits.BinTableHDU`_ object.

        Returns:
            `astropy.io.fits.BinTableHDU`_: HDU with the data.
        """
        cols = []
        hdr = fits.Header()
        for key in self.keys():
            if self[key] is None:
                continue
            # TODO: may be better to do this
            #   isinstance(self[key], (collections.Sequence, np.ndarray)):
            # This ignores the defined otype...
            if isinstance(self[key], (list, np.ndarray)):
                cols += [fits.Column(name=key,
                                     format=io.rec_to_fits_type(self[key], single_row=True),
                                     dim=io.rec_to_fits_col_dim(self[key], single_row=True),
                                     array=np.expand_dims(self[key], 0))]
            else:
                hdr[key.upper()] = (self[key], self.datamodel[key]['descr'])

        return fits.BinTableHDU.from_columns(cols, header=hdr, name=name)

    @classmethod
    def from_hdu(cls, hdu):
        """
        Instantiate the object from an `astropy.io.fits.BinTableHDU`_.

        Args:
            hdu (`astropy.io.fits.BinTableHDU`_):
                HDU containing the data to read.
        """
        d = dict.fromkeys(cls.datamodel.keys())
        for key in cls.datamodel.keys():
            if key.upper() in hdu.header.keys():
                d[key] = hdu.header[key.upper()]
                continue
            if key in hdu.columns.names:
                # NOTE: 0 is because all data are saved in a single row
                d[key] = hdu.data[key][0]
                continue

        # NOTE: We can't use `cls(d)` here because this will call the
        # `__init__` method of the derived class and we need to use the
        # `__init__` of the base class instead. So below, I get an
        # empty instance of the derived class using `__new__`, call the
        # parent `__init__`, and then return the result. This is the
        # first time I've used `__new__` so we may want to tread
        # carefully, but it seems to work.
        self = super().__new__(cls)
        super(cls, self).__init__(d)
        return self

    @classmethod
    def from_file(cls, ifile, ext):
        """
        Instantiate the object from an extension in the specified fits file.

        This is a convenience wrapper for :func:`from_hdu`.
        
        Args:
            ifile (:obj:`str`):
                Fits file with the data to read
            ext (:obj:`str`, :obj:`int`):
                Name or index of the extension in the fits file with
                the data.

        Raises:
            FileNotFoundError:
                Raised if the specified file does not exist.
        """
        if not os.path.isfile(ifile):
            raise FileNotFoundError('{0} does not exist!'.format(ifile))
        with fits.open(ifile) as hdu:
            return cls.from_hdu(hdu[ext])

