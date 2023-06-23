"""
Implements classes and function for the PypeIt data model.

.. _data-container:

DataContainer
-------------

:class:`DataContainer` objects provide a utility for
enforcing a specific datamodel on an object, and provides convenience
routines for writing data to fits files. The class itself is an
abstract base class that cannot be directly instantiated. As a base
class, :class:`DataContainer` objects are versatile, but they have
their limitations.

Derived classes must do the following:

    - Define a class attribute called ``datamodel``. See the examples
      below for their format.
    - Provide an :func:`__init__` method that defines the
      instantiation calling sequence and passes the relevant
      dictionary to this base-class instantiation.
    - Provide a :func:`_validate` method, if necessary, that
      processes the data provided in the `__init__` into a complete
      instantiation of the object. This method and the
      :func:`__init__` method are the *only* places where attributes
      can be added to the class.
    - Provide a :func:`_bundle` method that reorganizes the datamodel
      into partitions that can be written to one or more fits
      extensions. More details are provided in the description of
      :func:`DataContainer._bundle`.
    - Provide a :func:`_parse` method that parses information in one
      or more fits extensions into the appropriate datamodel. More
      details are provided in the description of
      :func:`DataContainer._parse`.

.. note::

    The attributes of the class are *not required* to be a part of
    the ``datamodel``; however, it makes the object simpler when they
    are. Any attributes that are not part of the ``datamodel`` must
    be defined in either the :func:`__init__` or :func:`_validate`
    methods; otherwise, the class with throw an ``AttributeError``.

Here are some examples of how to and how not to use them.

Defining the datamodel
++++++++++++++++++++++

The format of the ``datamodel`` needed for each implementation of a
:class:`DataContainer` derived class is as follows.

The datamodel itself is a class attribute (i.e., it is a member of
the class, not just of an instance of the class). The datamodel is a
dictionary of dictionaries: Each key of the datamodel dictionary
provides the name of a given datamodel element, and the associated
item (dictionary) for the datamodel element provides the type and
description information for that datamodel element. For each
datamodel element, the dictionary item must provide:

    - ``otype``: This is the type of the object for this datamodel item. E.g.,
      for a float or a `numpy.ndarray`_, you would set ``otype=float`` and
      ``otype=np.ndarray``, respectively.  The ``otype`` can also be a tuple
      with optional types.  Beware optional types that are themselves
      DataContainers.  This works for
      :class:`~pypeit.images.pypeitimage.PypeItImage`, but it likely needs more
      testing.

    - ``descr``: This provides a text description of the datamodel
      element. This is used to construct the datamodel tables in the
      pypeit documentation.

If the object type is a `numpy.ndarray`_, you should also provide the
``atype`` keyword that sets the type of the data contained within the
array. E.g., for a floating point array containing an image, your
datamodel could be simply::

    datamodel = {'image' : dict(otype=np.ndarray, atype=float, descr='My image')}

More advanced examples are given below.

Basic container
+++++++++++++++

Here's how to create a derived class for a basic container that
holds two arrays and a metadata parameter::

    import numpy as np
    import inspect

    from pypeit.datamodel import DataContainer

    class BasicContainer(DataContainer):
        datamodel = {'vec1': dict(otype=np.ndarray, atype=float, descr='Test'),
                     'meta1': dict(otype=str, decr='test'),
                     'arr1': dict(otype=np.ndarray, atype=float, descr='test')}

        def __init__(self, vec1, meta1, arr1):
            # All arguments are passed directly to the container
            # instantiation
            args, _, _, values = inspect.getargvalues(inspect.currentframe())
            super(BasicContainer, self).__init__({k: values[k] for k in args[1:]}) 

        def _bundle(self):
            # Use the base class _bundle Specify the extension
            return super(BasicContainer, self)._bundle(ext='basic')

With this implementation:

    - You can instantiate the container so that number of data table
      rows would be the same (10)::

        data = BasicContainer(np.arange(10), 'length=10', np.arange(30).reshape(10,3))

    - After instantiating, access the data like this::

        # Get the data model keys
        keys = list(data.keys())

        # Access the datamodel as attributes ...
        print(data.vec1)
        # ... or as items
        print(data['meta1'])

    - Attributes and items are case-sensitive::

        # Faults because of case-sensitive attributes/items
        print(data.Vec1)
        print(data['MeTa1'])

    - The attributes of the container can only be part of the
      datamodel or added in either the :func:`DataContainer.__init__`
      or :func:`DataContainer._validate` methods::

        # Faults with KeyError
        data['newvec'] = np.arange(10)
        test = data['newvec']

        # Faults with AttributeError
        data.newvec = np.arange(10)
        test = data.newvec

    - The :class:`DataContainer` also enforces strict types for
      members of the datamodel::

        # Faults with TypeError
        data.vec1 = 3
        data.meta1 = 4.

    - Read/Write the data from/to a fits file. In this instantiation,
      the data is written to a `astropy.io.fits.BinTableHDU` object;
      the table has 10 rows because the shape of the arrays match
      this. The file I/O routines look like this::

        # Write to a file
        ofile = 'test.fits'
        data.to_file(ofile)

        # Write to a gzipped file
        ofile = 'test.fits.gz'
        data.to_file(ofile)

        # Test written data against input
        with fits.open(ofile) as hdu:
            print(len(hdu))
            # 2: The primary extension and a binary table with the data

            print(hdu[1].name)
            # BASIC: As set by the _bundle method
            
            print(len(hdu['BASIC'].data))
            # 10: The length of the data table
            
            print(hdu['BASIC'].columns.names)
            # ['vec1', 'arr1']: datamodel elements written to the table columns
            
            print(hdu['BASIC'].header['meta1'])
            # 'length=10': int, float, or string datamodel components are written to headers

    - If the shape of the first axis of the arrays (number of rows)
      do not match, the arrays are written as single elements of a
      table with one row::

        # Number of rows are mismatched 
        data = BasicContainer(np.arange(10), 'length=1', np.arange(30).reshape(3,10))
        data.to_file(ofile)

        with fits.open(ofile) as hdu:
            print(len(hdu))
            # 2: The primary extension and a binary table with the data
            print(len(hdu['BASIC'].data))
            # 1: All of the data is put into a single row

Mixed Object Containers
+++++++++++++++++++++++

:class:`DataContainer` objects can also contain multiple arrays
and/or `astropy.table.Table`_ objects. However, multiple Tables or
combinations of arrays and Tables cannot be bundled into individual
extensions. Here are two implementations of a mixed container, a good
one and a bad one::

    import numpy as np
    import inspect

    from astropy.table import Table

    from pypeit.datamodel import DataContainer

    class GoodMixedTypeContainer(DataContainer):
        datamodel = {'tab1': dict(otype=Table, descr='Test'),
                     'tab1len': dict(otype=int, descr='test'),
                     'arr1': dict(otype=np.ndarray, descr='test'),
                     'arr1shape': dict(otype=tuple, descr='test')}

        def __init__(self, tab1, arr1):
            # All arguments are passed directly to the container
            # instantiation, but the list is incomplete
            args, _, _, values = inspect.getargvalues(inspect.currentframe())
            super(GoodMixedTypeContainer, self).__init__({k: values[k] for k in args[1:]}) 

        def _validate(self):
            # Complete the instantiation
            self.tab1len = len(self.tab1)
            # NOTE: DataContainer does allow for tuples, but beware because
            # they have to saved to the fits files by converting them to strings
            # and writing them to the fits header.  So the tuples should be
            # short!  See _bundle below, and in DataContainer.
            self.arr1shape = self.arr1.shape

        def _bundle(self):
            # Bundle so there's only one Table and one Image per extension
            return [{'tab1len': self.tab1len, 'tab1': self.tab1},
                    {'arr1shape': str(self.arr1shape), 'arr1': self.arr1}]


    class BadMixedTypeContainer(GoodMixedTypeContainer):
        def _bundle(self):
            # Use default _bundle method, which will try to put both tables
            # in the same extension.  NOTE: Can't use super here because
            # GoodMixedTypeContainer doesn't have an 'ext' argument
            return DataContainer._bundle(self, ext='bad')

With this implementation:

    - To instantiate::

        x = np.arange(10)
        y = np.arange(10)+5
        z = np.arange(30).reshape(10,3)

        arr1 = np.full((3,3,3), -1)
        tab1 = Table(data=({'x':x,'y':y,'z':z}), meta={'test':'this'})

        data = GoodMixedTypeContainer(tab1, arr1)

    - Data access::

        print(data.tab1.keys())
        # ['x', 'y', 'z']
        print(assert data.tab1.meta['test'])
        # 'this'
        print(data.arr1shape)
        # (3,3,3)

    - Construct an `astropy.io.fits.HDUList`::

        hdu = data.to_hdu(add_primary=True)

        print(len(hdu))
        # 3: Includes the primary HDU, one with the table, and one with the array

        print([h.name for h in hdu])
        # ['PRIMARY', 'TAB1', 'ARR1']

    - Tuples are converted to strings::

        print(hdu['ARR1'].header['ARR1SHAPE'])
        # '(3,3,3)'

    - The tuples are converted back from strings when they're read
      from the HDU::

        _data = GoodMixedTypeContainer.from_hdu(hdu)
        print(_data.arr1shape)
        # (3,3,3)

    - Table metadata is also written to the header, which can be
      accessed with case-insensitive keys::

        print(hdu['TAB1'].header['TEST'])
        print(hdu['TAB1'].header['TesT'])
        # Both print: 'this'

    - However, it's important to note that the keyword case gets
      mangled when you read it back in. This has to do with the
      Table.read method and I'm not sure there's anything we can do
      about it without the help of the astropy folks. We recommend
      table metadata use keys that are in all caps::

        # Fails because of the string case
        print(_data.tab1.meta['test'])
        # This is okay
        print(_data.tab1.meta['TEST'])

    - The difference between the implementation of
      ``BadMixedTypeContainer`` and ``GoodMixedTypeContainer`` has to
      do with how the data is bundled into HDU extensions. The
      ``BadMixedTypeContainer`` will instantiate fine::

        data = BadMixedTypeContainer(tab1, arr1)
        print(data.tab1.keys())
        # ['x', 'y', 'z']
        print(data.tab1.meta['test'])
        # 'this'

    - But it will barf when you try to reformat/write the data
      because you can't write both a Table and an array to a single
      HDU::

        # Fails
        hdu = data.to_hdu()

Complex Instantiation Methods
+++++++++++++++++++++++++++++

All of the :class:`DataContainer` above have had simple instatiation
methods. :class:`DataContainer` can have more complex instantiation
methods, but there are significant limitations to keep in made.
Consider::

    class BadInitContainer(DataContainer):
        datamodel = {'inp1': dict(otype=np.ndarray, descr='Test'),
                     'inp2': dict(otype=np.ndarray, descr='test'),
                     'out': dict(otype=np.ndarray, descr='test'),
                     'alt': dict(otype=np.ndarray, descr='test')}

        def __init__(self, inp1, inp2, func='add'):
            args, _, _, values = inspect.getargvalues(inspect.currentframe())
            super(BadInitContainer, self).__init__({k: values[k] for k in args[1:]}) 


    class DubiousInitContainer(DataContainer):
        datamodel = {'inp1': dict(otype=np.ndarray, descr='Test'),
                     'inp2': dict(otype=np.ndarray, descr='test'),
                     'out': dict(otype=np.ndarray, descr='test'),
                     'alt': dict(otype=np.ndarray, descr='test')}

        def __init__(self, inp1, inp2, func='add'):
            # If any of the arguments of the init method aren't actually
            # part of the datamodel, you can't use the nominal two lines
            # used in all the other examples above.  You have to be specific
            # about what gets passed to super.__init__.  See the
            # BadInitContainer example above and the test below.  WARNING:
            # I'm not sure you would ever want to do this because it can
            # lead to I/O issues; see the _validate function.
            self.func = func
            super(DubiousInitContainer, self).__init__({'inp1': inp1, 'inp2':inp2})

        def _validate(self):
            # Because func isn't part of the data model, it won't be part of
            # self if the object is instantiated from a file.  So I have to
            # add it here.  But I don't know what the value of the attribute
            # was for the original object that was written to disk.  This is
            # why you likely always want anything that's critical to setting
            # up the object to be part of the datamodel so that it gets
            # written to disk.  See the testing examples for when this will
            # go haywire.
            if not hasattr(self, 'func'):
                self.func = None
            if self.func not in [None, 'add', 'sub']:
                raise ValueError('Function must be either \'add\' or \'sub\'.')

            # This is here because I don't want to overwrite something that
            # might have been read in from an HDU, particularly given that
            # func will be None if reading from a file!
            if self.out is None:
                print('Assigning out!')
                if self.func is None:
                    raise ValueError('Do not know how to construct out attribute!')
                self.out = self.inp1 + self.inp2 if self.func == 'add' else self.inp1 - self.inp2

        # I'm not going to overwrite _bundle, so that the nominal approach
        # is used.


    class ComplexInitContainer(DataContainer):
        datamodel = {'inp1': dict(otype=np.ndarray, descr='Test'),
                     'inp2': dict(otype=np.ndarray, descr='test'),
                     'out': dict(otype=np.ndarray, descr='test'),
                     'alt': dict(otype=np.ndarray, descr='test'),
                     'func': dict(otype=str, descr='test')}

        def __init__(self, inp1, inp2, func='add'):
            # Since func is part of the datamodel now, we can use the normal
            # two intantiation lines.
            args, _, _, values = inspect.getargvalues(inspect.currentframe())
            super(ComplexInitContainer, self).__init__({k: values[k] for k in args[1:]}) 

        def _validate(self):
            if self.func not in ['add', 'sub']:
                raise ValueError('Function must be either \'add\' or \'sub\'.')
            # This is here because I don't want to overwrite something that
            # might have been read in from an HDU, even though they should
            # nominally be the same!
            if self.out is None:
                print('Assigning out!')
                if self.func is None:
                    raise ValueError('Do not know how to construct out attribute!')
                self.out = self.inp1 + self.inp2 if self.func == 'add' else self.inp1 - self.inp2


With this implementation:

    - The instantiation of the ``BadInitContainer`` will fail because
      the init arguments all need to be part of the datamodel for it
      to work::

        x = np.arange(10)
        y = np.arange(10)+5

        data = BadInitContainer(x,y)
        # Fails with AttributeError

    - The following instantiation is fine because
      ``DubiousInitContainer`` handles the fact that some of the
      arguments to :func:`__init__` are not part of the datamodel::

        data = DubiousInitContainer(x,y)
        print(np.array_equal(data.out, data.inp1+data.inp2))
        # True

    - One component of the data model wasn't instantiated, so it will
      be None::
        
        print(data.alt is None)
        # True

    - The problem with ``DubiousInitContainer`` is that it has
      attributes that cannot be reinstantiated from what's written to
      disk. That is, :class:`DataContainers` aren't really fully
      formed objects unless all of its relevant attributes are
      components of the data model::

        _data = DubiousInitContainer.from_hdu(hdu)
        print(_data.func == DubiousInitContainer(x,y).func)
        # False
    
    - This is solved by adding func to the datamodel::
    
        data = ComplexInitContainer(x,y)
        _data = ComplexInitContainer.from_hdu(data.to_hdu(add_primary=True))
        print(data.func == _data.func)
        # True

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
from pathlib import Path

from IPython import embed

import numpy as np
import inspect

from astropy.io import fits
from astropy.table import Table

from pypeit import io
from pypeit import msgs

# TODO: There are methods in, e.g., doc/scripts/build_specobj_rst.py that output
# datamodels for specific datacontainers.  It would be useful if we had
# generalized methods in this base class that
#   - construct an rst table that documents the datamodel.
#   - construct an rst table that documents each Table component of the datamodel
#     (like a nested DataContainer)

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
          instantiation calling sequence and passes the relevant
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

    .. note::

        The attributes of the class are *not required* to be a part
        of the ``datamodel``; however, it makes the object simpler
        when they are. Any attributes that are not part of the
        ``datamodel`` must be defined in either the :func:`__init__`
        or :func:`_validate` methods; otherwise, the class with throw
        an ``AttributeError``.

    .. todo::

        Add a copy method

    Args:
        d (:obj:`dict`, optional):
            Dictionary to copy to the internal attribute dictionary.
            All of the keys in the dictionary *must* be elements of
            the ``datamodel``. Any attributes that are not part of
            the ``datamodel`` can be set in the :func:`__init__` or
            :func:`_validate` methods of a derived class. If None,
            the object is instantiated with all of the relevant data
            model attributes but with all of those attributes set to
            None.

    """

    # Define the class version
    version = None
    """
    Provides the string representation of the class version.

    This is currently put to minimal use so far, but will used for
    I/O verification in the future.

    Each derived class should provide a version to guard against data
    model changes during development.
    """

    hdu_prefix = None
    """
    If set, all HDUs generated for this DataContainer will have this
    prefix. This can be set independently for each DataContainer
    derived class; however, it always defaults to None for the base
    class. Be wary of nested DataContainer's!!
    """

    output_to_disk = None
    """
    If set, this limits the HDU extensions that are written to the
    output file. Note this is the name of the extension, including
    the hdu_prefix, not necessarily the names of specific datamodel
    components.
    """

    one_row_table = False
    """
    Force the full datamodel to be encapsulated into an `astropy.table.Table`_
    with a single row when written to disk.  Beware that this requires that this
    is possible!  See, e.g.,
    :class:`~pypeit.images.detector_container.DetectorContainer`.
    """

    # Define the data model
    datamodel = None
    """
    Provides the class data model. This is undefined in the abstract
    class and should be overwritten in the derived classes.

    The format of the ``datamodel`` needed for each implementation of
    a :class:`DataContainer` derived class is as follows.

    The datamodel itself is a class attribute (i.e., it is a member
    of the class, not just of an instance of the class). The
    datamodel is a dictionary of dictionaries: Each key of the
    datamodel dictionary provides the name of a given datamodel
    element, and the associated item (dictionary) for the datamodel
    element provides the type and description information for that
    datamodel element. For each datamodel element, the dictionary
    item must provide:

        - ``otype``: This is the type of the object for this
          datamodel item. E.g., for a float or a `numpy.ndarray`_,
          you would set ``otype=float`` and ``otype=np.ndarray``,
          respectively.

        - ``descr``: This provides a text description of the
          datamodel element. This is used to construct the datamodel
          tables in the pypeit documentation.

    If the object type is a `numpy.ndarray`_, you should also provide
    the ``atype`` keyword that sets the type of the data contained
    within the array. E.g., for a floating point array containing an
    image, your datamodel could be simply::

        datamodel = {'image' : dict(otype=np.ndarray, atype=float, descr='My image')}

    More advanced examples are given in the top-level module documentation.

    Currently, ``datamodel`` components are restricted to have
    ``otype`` that are :obj:`tuple`, :obj:`int`, :obj:`float`,
    ``numpy.integer``, ``numpy.floating``, `numpy.ndarray`_, or
    `astropy.table.Table`_ objects. E.g., ``datamodel`` values for
    ``otype`` *cannot* be :obj:`dict`.
    """
    # TODO: Enable multiple possible types for the datamodel elements?
    # I.e., allow `otype` to be a tuple of the allowed object types? It
    # looks like this is already possible at least for some types, see
    # pypeit.tracepca.TracePCA.reference_row.

    internals = None
    """
    A list of strings identifying a set of internal attributes that are not part
    of the datamodel.
    """

    def __init__(self, d=None):

        # Data model must be defined
        if self.datamodel is None:
            raise ValueError('Data model for {0} is undefined!'.format(self.__class__.__name__))
        if self.version is None:
            raise ValueError('Must define a version for the class.')

        # Ensure the dictionary has all the expected keys
        self.__dict__.update(dict.fromkeys(self.datamodel.keys()))

        # Initialize other internals
        self._init_internals()

        # Finalize the instantiation.
        # NOTE: The key added to `__dict__` by this call is always
        # `_DataContainer__initialised`, regardless of whether or not
        # the call to this `__init__` is from the derived class. This
        # is why I can check for `_DataContainer__initialised` is in
        # `self.__dict__`, even for derived classes. But is there a way
        # we could just check a boolean instead?
        self._init_key = '_DataContainer__initialised'
        self.__initialised = True

        # Include the provided data and build-out the data model, if
        # data were provided
        if d is not None:

            # Input dictionary cannot have keys that do not exist in
            # the data model
            if not np.all(np.isin(list(d.keys()), list(self.datamodel.keys()))):
                raise AttributeError('Coding error: Initialization arguments do not match '
                                     'data model!')

            # Assign the values provided by the input dictionary
            #self.__dict__.update(d)  # This by-passes the data model checking

            ## Assign the values provided by the input dictionary
            for key in d.keys():

                # DataContainer element can be one of multiple DataContainer types
                optional_dc_type = isinstance(self.datamodel[key]['otype'], tuple) \
                        and any([obj_is_data_container(t) for t in self.datamodel[key]['otype']])
                if optional_dc_type and isinstance(d[key], dict):
                    for t in self.datamodel[key]['otype']:
                        try:
                            dc = t(**d[key])
                        except:
                            dc = None
                            continue
                        else:
                            break
                    if dc is None:
                        msgs.error(f'Could not assign dictionary element {key} to datamodel '
                                   f'for {self.__class__.__name__}.', cls='PypeItDataModelError')
                    setattr(self, key, dc)
                    continue
                
                # Nested DataContainer?
                if obj_is_data_container(self.datamodel[key]['otype']) and isinstance(d[key], dict):
                    setattr(self, key, self.datamodel[key]['otype'](**d[key]))
                    continue

                setattr(self, key, d[key])

        # Validate the object
        self._validate()

    def _init_internals(self):
        """
        Add internal variables to the object before initialization completes

        These should be set to None
        """
        if self.internals is None:
            return
        for attr in self.internals:
            setattr(self, attr, None) 

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
              :func:`~pypeit.io.dict_to_hdu` for dictionaries with
              single array or `astropy.table.Table`_ items).

        The item to be written can be a single array for an
        `astropy.io.fits.ImageHDU`_, an `astropy.table.Table`_ for a
        `astropy.io.fits.BinTableHDU`_, or a dictionary; see
        :func:`~pypeit.io.write_to_hdu`.

        For how these objects are parsed into the HDUs, see
        :func:`~pypeit.datamodel.DataContainer.to_hdu`.

        The default behavior implemented by this base class just parses the
        attributes into a single dictionary based on the datamodel, under the
        expectation that it is all written to a single fits extension. Note that
        this will fault if the datamodel contains:

            - a dictionary object
            - more than one `astropy.table.Table`_, or
            - the combination of an `astropy.table.Table`_ and one or more
              array-like objects (:obj:`list` or `numpy.ndarray`_)

        Certain **restrictions** apply to how the data can be bundled
        for the general parser implementation (:func:`_parse`) to
        work correctly. These restrictions are:

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
        """
        d = {}
        for key in self.keys():
            if self[key] is not None and transpose_arrays \
                    and self.datamodel[key]['otype'] == np.ndarray:
                d[key] = self[key].T
            elif self.datamodel[key]['otype'] == tuple:
                # TODO: Anything with tuple type that is None will be
                # converted to 'None'. Is that what we want, or do we
                # want to set it to None so that it's not written?
                d[key] = str(self[key])
            else:
                d[key] = self[key]

        if self.one_row_table:
            # Attempt to stuff the full datamodel into an astropy Table with a
            # single row
            # NOTE: This is annoyingly complicated to avoid adding elements that
            # are None to the table.
            del_keys = []
            for key, val in d.items():
                if isinstance(val, np.ndarray):
                    d[key] = np.expand_dims(val, 0) 
                elif val is None:
                    del_keys += [key]
                else:
                    d[key] = [val]
            for key in del_keys:
                del d[key]

            try:
                d = Table(d)
            except:
                msgs.error(f'Cannot force all elements of {self.__class__.__name__} datamodel'
                           'into a single-row astropy Table!', cls='PypeItDataModelError')

        return [d] if ext is None else [{ext:d}]

    @classmethod
    def confirm_class(cls, name, allow_subclasses=False):
        """
        Confirm that the provided name of a datamodel class matches the calling
        class.

        Args:
            name (:obj:`str`):
                The name of the class to check against.
            allow_subclasses (:obj:`bool`, optional):
                Allow subclasses of the calling class to be included in the
                check.
        Returns:
            :obj:`bool`: Class match confirmed if True
        """
        is_exact = name == cls.__name__
        if is_exact or not allow_subclasses:
            return is_exact
        # NOTE: Unlike pypeit.utils.all_subclasses, this does *not* recursively
        # check for all the subclasses.  We can do that, but it's not necessary
        # yet.
        return is_exact or name in [c.__name__ for c in cls.__subclasses__()]

    @classmethod
    def _parse(cls, hdu, ext=None, ext_pseudo=None, transpose_table_arrays=False, 
               hdu_prefix=None, allow_subclasses=False):
        """
        Parse data read from one or more HDUs.

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
              also be included in the HDUs constructed by
              :func:`to_hdu`. However, the
              `astropy.table.Table.read`_ method always returns the
              metadata with capitalized keys. This means that,
              regardless of the capitalization of the metadata
              keywords when the data is written, **they will be
              upper-case when read by this function**!

            - Use ``allow_subclasses`` with care!  The expectation is that all
              the subclasses have exactly the same datamodel and version number
              as the main class.

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
            hdu (`astropy.io.fits.HDUList`_, `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_):
                The HDU(s) to parse into the instantiation dictionary.
            ext (:obj:`int`, :obj:`str`, :obj:`list`, optional):
                One or more extensions with the data. If None, the
                function trolls through the HDU(s) and parses the
                data into the datamodel.
            ext_pseudo (:obj:`int`, :obj:`str`, :obj:`list`, optional):
                Pseudonyms to use for the names of the HDU extensions instead of
                the existing extension names.  Similar to the use of
                ``hdu_prefix``, this is useful for parsing nested
                :class:`DataContainer` objects; i.e., a child
                :class:`DataContainer` may be assigned to an extension based on
                the datamodel of its parent, meaning that it can't be properly
                parsed by the parent's method.  If used, the number of extension
                pseudonyms provided **must** match the list returned by::
                    
                    io.hdu_iter_by_ext(hdu, ext=ext, hdu_prefix=hdu_prefix)

            transpose_table_arrays (:obj:`bool`, optional):
                Tranpose *all* the arrays read from any binary
                tables. This is meant to invert the use of
                ``transpose_arrays`` in :func:`_bound`.
            hdu_prefix (:obj:`str`, optional):
                Only parse HDUs with extension names matched to this prefix. If
                None, :attr:`ext` is used. If the latter is also None, all HDUs
                are parsed. See :func:`pypeit.io.hdu_iter_by_ext`.
            allow_subclasses (:obj:`bool`, optional):
                Allow subclasses of the calling class to successfully pass the
                check on whether or not it can parse the provided HDUs.  That
                is, if the class provided in the HDU header is "A", setting this
                parameter to True means that parsing will continue as long as
                "A" is a subclass of the calling class.  In ``PypeIt`` this is
                needed, e.g., for :class:`~pypeit.sensfunc.SensFunc` to
                successfully parse files written by either of its subclasses,
                :class:`~pypeit.sensfunc.IRSensFunc` or 
                :class:`~pypeit.sensfunc.UVISSensFunc`.  This is possible
                because the datamodel is defined by the parent class and not
                altered by the subclasses!  **Use with care! See function
                warnings.**
                
        Returns:
            :obj:`tuple`: Return four objects: (1) a dictionary with the
            datamodel contents, (2) a boolean flagging if the datamodel version
            checking has passed, (3) a boolean flagging if the datamodel type
            checking has passed, and (4) the list of parsed HDUs.

        """
        dm_version_passed = True
        dm_type_passed = True

        # Setup to iterate through the provided HDUs
        _ext, _hdu = io.hdu_iter_by_ext(hdu, ext=ext, hdu_prefix=hdu_prefix)
        _ext = np.atleast_1d(np.array(_ext, dtype=object))

        _ext_pseudo = _ext if ext_pseudo is None else np.atleast_1d(ext_pseudo)

        if len(_ext_pseudo) != len(_ext):
            msgs.error(f'Length of provided extension pseudonym list must match number of '
                       f'extensions selected: {len(_ext)}.', cls='PypeItDataModelError')

        str_ext = np.logical_not([isinstance(e, (int, np.integer)) for e in _ext_pseudo])

        # Construct instantiation dictionary
        _d = dict.fromkeys(cls.datamodel.keys())

        # Log if relevant data is found for this datamodel, allowing for
        # DataContainers that have no data, although such a usage case should be
        # rare.
        if np.all([_hdu[e].data is None for e in _ext]):
            msgs.warn(f'Extensions to be read by {cls.__name__} have no data!')
            # This is so that the returned booleans for reading the
            # data are not tripped as false!
            found_data = True
        else:
            found_data = False

        # NOTE: The extension and keyword comparisons are complicated
        # because the fits standard is to force these all to be
        # capitalized, while the datamodel doesn't
        # implement this restriction.

        # Handle hdu_prefix
        if hdu_prefix is not None:
            prefix = hdu_prefix
        else:
            prefix = '' if cls.hdu_prefix is None else cls.hdu_prefix

        # Save the list of hdus that have been parsed
        parsed_hdus = []

        # HDUs can have dictionary elements directly.
        keys = np.array(list(_d.keys()))
        prefkeys = np.array([prefix+key.upper() for key in keys])
        indx = np.isin(prefkeys, _ext_pseudo[str_ext])

        if np.any(indx):
            found_data = True
            for e in keys[indx]:
                pseudo_hduindx = prefix+e.upper()
                hduindx = _ext[np.where(_ext_pseudo == pseudo_hduindx)[0][0]]
                # Add it to the list of parsed HDUs
                parsed_hdus += [hduindx]

                # Parse an extension that can be one of many DataContainers
                if isinstance(cls.datamodel[e]['otype'], tuple) \
                        and any([obj_is_data_container(t) for t in cls.datamodel[e]['otype']]):
                    for t in cls.datamodel[e]['otype']:
                        # TODO: May need a try block here...
                        data, dvpassed, dtpassed, _ = t._parse(_hdu[hduindx])
                        if dvpassed and dtpassed:
                            break
                    # Save the successfully parsed data, or the last one that
                    # failed and continue
                    _d[e] = data
                    dm_version_passed &= dvpassed
                    dm_type_passed &= dtpassed
                    continue

                # Parse an extension that can only be one DataContainer
                if obj_is_data_container(cls.datamodel[e]['otype']):
                    # Parse the DataContainer
                    # TODO: This only works with single extension
                    # DataContainers. Do we want this to be from_hdu
                    # instead and add chk_version to _parse?
                    # TODO: Should this propagate allow_subclasses?  It does
                    # not, for now.
                    _d[e], p1, p2, _ = cls.datamodel[e]['otype']._parse(_hdu[hduindx])
                    dm_version_passed &= p1
                    dm_type_passed &= p2
                    continue

                # Parse Image or Table data
                dm_type_passed &= cls.confirm_class(_hdu[hduindx].header['DMODCLS'],
                                                    allow_subclasses=allow_subclasses)
                dm_version_passed &= _hdu[hduindx].header['DMODVER'] == cls.version
                # Grab it
                _d[e] = _hdu[hduindx].data if isinstance(_hdu[hduindx], fits.ImageHDU) \
                    else Table.read(_hdu[hduindx]).copy()

        for e in _ext:
            if 'DMODCLS' not in _hdu[e].header.keys() or 'DMODVER' not in _hdu[e].header.keys() \
                    or not cls.confirm_class(_hdu[e].header['DMODCLS'],
                                             allow_subclasses=allow_subclasses):
                # Can't be parsed
                continue
            # Check for header elements, but do not over-ride existing items
            indx = np.isin([key.upper() for key in keys], list(_hdu[e].header.keys()))
            if np.any(indx):
                found_data = True
                parsed_hdus += [e if _hdu[e].name == '' else _hdu[e].name]
                dm_type_passed &= cls.confirm_class(_hdu[e].header['DMODCLS'],
                                                    allow_subclasses=allow_subclasses)
                dm_version_passed &= _hdu[e].header['DMODVER'] == cls.version
                for key in keys[indx]:
                    if key in _d.keys() and _d[key] is not None:
                        continue
                    _d[key] = _hdu[e].header[key.upper()] if cls.datamodel[key]['otype'] != tuple \
                                else eval(_hdu[e].header[key.upper()])
            if isinstance(e, (str, np.str_)) and e in prefkeys:
                # Already parsed this above
                continue
            # Parse BinTableHDUs
            if isinstance(_hdu[e], fits.BinTableHDU) \
                    and np.any(np.isin(list(cls.datamodel.keys()), _hdu[e].columns.names)):
                parsed_hdus += [e if _hdu[e].name is None else _hdu[e].name]
                found_data = True
                # Datamodel checking
                dm_type_passed &= cls.confirm_class(_hdu[e].header['DMODCLS'],
                                                    allow_subclasses=allow_subclasses)
                dm_version_passed &= _hdu[e].header['DMODVER'] == cls.version
                # If the length of the table is 1, assume the table
                # data had to be saved as a single row because of shape
                # differences.
                single_row = len(_hdu[e].data) == 1
                for key in _hdu[e].columns.names:
                    if key in cls.datamodel.keys():
                        _d[key] = _hdu[e].data[key][0] \
                                        if (single_row and _hdu[e].data[key].ndim > 1) \
                                        else _hdu[e].data[key]
                        if transpose_table_arrays:
                            _d[key] = _d[key].T

        # Two annoying hacks:
        #   - Hack to expunge charray which are basically deprecated and
        #     cause trouble.
        #   - Hack to force native byte ordering
        for key in _d:
            if isinstance(_d[key], np.chararray):
                _d[key] = np.asarray(_d[key])
            elif isinstance(_d[key], np.ndarray) and _d[key].dtype.byteorder not in ['=', '|']:
                _d[key] = _d[key].astype(_d[key].dtype.type)

        # Unpack if necessary
        if cls.one_row_table:
            # WARNING: This has only been tested for single-extension
            # DataContainers!!

            # NOTE: This works for the 1D vector elements of DetectorContainer,
            # but it's untested for any higher dimensional arrays...
            _d = {key:val if cls.datamodel[key]['otype'] == np.ndarray or val is None else val[0] 
                    for key, val in _d.items()}
            # NOTE: Annoyingly, casting becomes awkward when forcing all the
            # data into an astropy Table.  E.g., int becomes np.int64, which
            # won't be accepted by the strict typing rules used when setting
            # attributes/items.  The following is a work-around for this
            # specific case, but it likely points to a larger issue that we may
            # need to address in DataContainer.
            for key in _d.keys():
                if _d[key] is None or key not in cls.datamodel:
                    continue
                types = cls.datamodel[key]['otype'] \
                            if isinstance(cls.datamodel[key]['otype'], (list, tuple)) \
                            else [cls.datamodel[key]['otype']]
                if any([isinstance(_d[key], t) for t in types]):
                    continue
                for t in types:
                    # NOTE: 'equiv' is too strict for windows?  Causes CI tests
                    # for windows to fail.
#                    if np.can_cast(_d[key], t, casting='equiv'):
                    if np.can_cast(_d[key], t, casting='same_kind'):
                        _d[key] = t(_d[key])
                        break

        # Return
        return _d, dm_version_passed and found_data, dm_type_passed and found_data, \
                    np.unique(parsed_hdus).tolist()


    def __getattr__(self, item):
        """Maps values to attributes.
        Only called if there *isn't* an attribute with this name
        """
        try:
            return self.__getitem__(item)
        except KeyError as e:
            raise AttributeError('{0} is not an attribute of {1}!'.format(item,
                                    self.__class__.__name__)) from e

    def __setattr__(self, item, value):
        """
        Set the attribute value.

        Attributes are restricted to those defined by the datamodel.
        """
        # version is immutable
        if item == 'version':
            raise TypeError('Internal version does not support assignment.')
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
                raise AttributeError('{0} is not part of the internals nor data model!'.format(item)) from e

    def __setitem__(self, item, value):
        """
        Access and set an attribute identically to a dictionary item.

        Items are restricted to those defined by the datamodel.
        """
        if item not in self.__dict__.keys():
            raise KeyError('Key {0} not part of the internals nor data model'.format(item))
        # Internal?
        if item not in self.keys():
            self.__dict__[item] = value
            return
        # Set datamodel item to None?
        if value is None:
            self.__dict__[item] = value
            return
        # Convert Path objects to string for saving in the datamodel
        if isinstance(value, Path):
            value = str(value)
        # Check data type
        if not isinstance(value, self.datamodel[item]['otype']):
            raise TypeError(f'Cannot assign object of type {type(value)} to {item}.\n'
                            f"Allowed type(s) are: {self.datamodel[item]['otype']}")
        # Array?
        if 'atype' in self.datamodel[item].keys():
            if not isinstance(value.flat[0], self.datamodel[item]['atype']):
                raise TypeError(f'Cannot assign array with data type {type(value.flat[0])} to '
                                f'{item} array.\nAllowed type(s) for the array are: '
                                f"{self.datamodel[item]['atype']}")
        # Set
        self.__dict__[item] = value

    def __getitem__(self, item):
        """Get an item directly from the internal dict."""
        try:
            return self.__dict__[item]
        except KeyError as e:
            raise KeyError(f'{item} is not an item in {self.__class__.__name__}.') from e

    def keys(self):
        """
        Return the keys for the data objects only

        Returns:
            :obj:`dict_keys`: The iterable with the data model keys.
        """
        return self.datamodel.keys()

    def check_populated(self, dm_items):
        """
        Check that a set of datamodel items are populated.

        Args:
            dm_items (:obj:`list`, :obj:`str`):
                One or more items in the datamodel to check.

        Returns:
            :obj:`bool`: Flag that *all* the requested datamodel items are
            populated (not None).
        """
        return np.all([key in self.keys() and self[key] is not None for key in dm_items])

    def _primary_header(self, hdr=None):
        """
        Construct a primary header that is included with the primary
        HDU extension produced by :func:`to_hdu`.

        Additional data can be added to the header for individual HDU
        extensions in :func:`to_hdu` as desired, but this function
        **should not** add any elements from the datamodel to the
        header.

        Args:
            hdr (`astropy.io.fits.Header`_, optional):
                Header for the primary extension. If None, set by
                :func:`pypeit.io.initialize_header()`.

        Returns:
            `astropy.io.fits.Header`_: Header object to include in
            the primary HDU.
        """
        return io.initialize_header() if hdr is None else hdr.copy()

    def _base_header(self, hdr=None):
        """
        Construct a base header that is included with all HDU
        extensions produced by :func:`to_hdu` (unless they are
        overwritten by a nested :class:`DataContainer`).

        Additional data can be added to the header for individual HDU
        extensions in :func:`to_hdu` as desired, but this function
        **should not** add any elements from the datamodel to the
        header.

        Args:
            hdr (`astropy.io.fits.Header`_, optional):
                Baseline header to add to all returned HDUs. If None,
                set by :func:`pypeit.io.initialize_header()`.

        Returns:
            `astropy.io.fits.Header`_: Header object to include in
            all HDU extensions.
        """
        # Copy primary header to all subsequent headers
        _hdr = self._primary_header(hdr=hdr)
        # Add DataContainer class name and datamodel version number.
        # This is not added to the primary header because single output
        # files can contain multiple DataContainer objects.
        _hdr['DMODCLS'] = (self.__class__.__name__, 'Datamodel class')
        _hdr['DMODVER'] = (self.version, 'Datamodel version')
        return _hdr

    @staticmethod
    def valid_write_to_hdu_type(obj):
        """
        Check if the provided object can be written to an `astropy.io.fits.HDUList`_.

        This needs to be consistent with :func:`~pypeit.io.write_to_hdu`.

        Args:
            obj (object):
                Object to write

        Returns:
            :obj:`bool`: Flag that object can be written.
        """
        return isinstance(obj, (DataContainer, dict, Table, np.ndarray, list))

    # TODO: Always have this return an HDUList instead of either that
    # or a normal list?
    # NOTE: This function should *not* include **kwargs.
    def to_hdu(self, hdr=None, add_primary=False, primary_hdr=None,
               force_to_bintbl=False, hdu_prefix=None):
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
            hdr (`astropy.io.fits.Header`, optional):
                Baseline header to add to all returned HDUs. If None,
                set by :func:`pypeit.io.initialize_header()`.
            add_primary (:obj:`bool`, optional):
                If False, the returned object is a simple
                :obj:`list`, with a list of HDU objects (either
                `astropy.io.fits.ImageHDU`_ or
                `astropy.io.fits.BinTableHDU`_). If true, the method
                constructs an `astropy.io.fits.HDUList` with a
                primary HDU, such that this call::

                    hdr = io.initialize_header()
                    hdu = fits.HDUList([fits.PrimaryHDU(header=primary_hdr)] + self.to_hdu(hdr=hdr))

                and this call::

                    hdu = self.to_hdu(add_primary=True)

                are identical.
            primary_hdr (`astropy.io.fits.Header`, optional):
                Header to add to the primary if add_primary=True
            force_to_bintbl (:obj:`bool`, optional):
                Force construction of a `astropy.io.fits.BinTableHDU`_ instead
                of an `astropy.io.fits.ImageHDU`_ when either there are no
                arrays or tables to write or only a single array is provided (as
                needed for, e.g., :class:`~pypeit.specobj.SpecObj`).  See
                :func:`~pypeit.io.write_to_hdu`.  
            hdu_prefix (:obj:`str`, optional):
                Prefix for the HDU names. If None, will use
                :attr:`hdu_prefix`. If the latter is also None, no
                prefix is added.

        Returns:
            :obj:`list`, `astropy.io.fits.HDUList`_: A list of HDUs, where the
            type depends on the value of ``add_primary``: If True, an
            `astropy.io.fits.HDUList`_ is returned, otherwise a :obj:`list` is
            returned.
        """
        # Bundle the data
        data = self._bundle()

        # Initialize the primary header (only used if add_primary=True)
        _primary_hdr = self._primary_header(hdr=primary_hdr)
        # Check that the keywords in the primary header do not overlap
        # with any datamodel keys.
        if _primary_hdr is not None:
            hdr_keys = np.array([k.upper() for k in self.keys()])
            indx = np.isin(hdr_keys, list(_primary_hdr.keys()))
            if np.sum(indx) > 1:
                msgs.error('CODING ERROR: Primary header should not contain keywords that are the '
                           'same as the datamodel for {0}.'.format(self.__class__.__name__),
                           cls='PypeItDataModelError')

        # Initialize the base header
        _hdr = self._base_header(hdr=hdr)
        # Check that the keywords in the base header do not overlap
        # with any datamodel keys.
        if _hdr is not None \
                and np.any(np.isin([k.upper() for k in self.keys()], list(_hdr.keys()))):
            msgs.error('CODING ERROR: Baseline header should not contain keywords that are the '
                       'same as the datamodel for {0}.'.format(self.__class__.__name__),
                       cls='PypeItDataModelError')

        # Construct the list of HDUs
        hdu = []
        for d in data:
            if isinstance(d, dict) and len(d) == 1 \
                    and DataContainer.valid_write_to_hdu_type(d[list(d.keys())[0]]):
                ext = list(d.keys())[0]
                # Allow for embedded DataContainer's
                if isinstance(d[ext], DataContainer):
                    _hdu = d[ext].to_hdu(add_primary=False)
                    # NOTE: The lines below allow extension names to be
                    # overridden by the input dictionary keyword for
                    # DataContainers that are confined to a single HDU.
                    # This is necessary because `hdu_prefix` must be a
                    # class attribute to allow for its inclusion in the
                    # read methods (e.g., `_parse`)
                    if len(_hdu) == 1 and _hdu[0].name != ext:
                        _hdu[0].name = ext
                    hdu += _hdu
                else:
                    hdu += [io.write_to_hdu(d[ext], name=ext, hdr=_hdr,
                                            force_to_bintbl=force_to_bintbl)]
            else:
                hdu += [io.write_to_hdu(d, hdr=_hdr, force_to_bintbl=force_to_bintbl)]

        # Prefixes
        _hdu_prefix = self.hdu_prefix if hdu_prefix is None else hdu_prefix
        if _hdu_prefix is not None:
            for ihdu in hdu:
                ihdu.name = f'{_hdu_prefix}{ihdu.name}'
        # Limit?
        if self.output_to_disk is not None:
            hdu = [h for h in hdu if h.name in self.output_to_disk]

        # Return
        return fits.HDUList([fits.PrimaryHDU(header=_primary_hdr)] + hdu) if add_primary else hdu

    @classmethod
    def from_hdu(cls, hdu, chk_version=True, **kwargs):
        """
        Instantiate the object from an HDU extension.

        This is primarily a wrapper for :func:`_parse`.

        Args:
            hdu (`astropy.io.fits.HDUList`_, `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_):
                The HDU(s) with the data to use for instantiation.
            chk_version (:obj:`bool`, optional):
                If True, raise an error if the datamodel version or
                type check failed. If False, throw a warning only.
            **kwargs:
                Passed directly to :func:`_parse`.
        """
        d, dm_version_passed, dm_type_passed, parsed_hdus = cls._parse(hdu, **kwargs)
        # Check version and type?
        if not dm_type_passed:
            msgs.error(f'The HDU(s) cannot be parsed by a {cls.__name__} object!',
                       cls='PypeItDataModelError')
        if not dm_version_passed:
            msg = f'Current version of {cls.__name__} object in code ({cls.version}) ' \
                  'does not match version used to write your HDU(s)!'
            if chk_version:
                msgs.error(msg, cls='PypeItDataModelError')
            else:
                msgs.warn(msg)

        # Instantiate
        # NOTE: We can't use `cls(d)`, where `d` is the dictionary returned by
        # `cls._parse`, because this will call the `__init__` method of the
        # derived class and we need to use the `__init__` of the base class
        # instead.  Instead, `d` is passed to `cls.from_dict`, which initiates
        # everything correctly.
        return cls.from_dict(d=d)

    @classmethod
    def from_dict(cls, d=None):
        """
        Instantiate from a dictionary.

        This is primarily to allow for instantiating classes from a
        file where the data has already been parsed. E.g., see how
        the :class:`~pypeit.tracepca.TracePCA` objects are
        instantiated in
        :class:`pypeit.edgetrace.EdgeTraceSet.from_hdu`. However,
        note that this does the bare minimum to instantiate the
        object. Any class-specific operations that are needed to
        complete the instantiation should be done by ovewritting this
        method; e.g., see :func:`pypeit.tracepca.TracePCA.from_dict`.

        Args:
            d (:obj:`dict`, optional):
                Dictionary with the data to use for instantiation.
        """
        # NOTE: We can't use `cls(d)` here because this will call the `__init__`
        # method of the derived class and we need to use the `__init__` of the
        # base class instead. So instead, I get an empty instance of the derived
        # class using `__new__`, call the parent `__init__`, and then return the
        # result. The call to `DataContainer.__init__` is explicit to deal with
        # objects inheriting from both DataContainer and other base classes.
        self = super().__new__(cls)
        DataContainer.__init__(self, d=d)
        return self

    def to_file(self, ofile, overwrite=False, checksum=True, **kwargs):
        """
        Write the data to a file.

        This is a convenience wrapper for :func:`to_hdu` and
        :func:`pypeit.io.write_to_fits`.  The ``add_primary`` parameter of
        :func:`to_hdu` is *always* true such that the first extension of the
        written fits file is *always* an empty primary header.

        Args:
            ofile (:obj:`str`, `Path`_):
                Fits file for the data. File names with '.gz'
                extensions will be gzipped; see
                :func:`pypeit.io.write_to_fits`.
            overwrite (:obj:`bool`, optional):
                Flag to overwrite any existing file.
            checksum (:obj:`bool`, optional):
                Passed to `astropy.io.fits.HDUList.writeto`_ to add
                the DATASUM and CHECKSUM keywords fits header(s).
            kwargs (:obj:`dict`, optional):
                Passed directly to :func:`to_hdu`.
        """
        # NOTE: This call does *not* need to also pass hdr to io.write_to_fits
        # because the first argument of the function is always an
        # astropy.io.fits.HDUList.
        io.write_to_fits(self.to_hdu(add_primary=True, **kwargs),
                         str(ofile), overwrite=overwrite, checksum=checksum)

    @classmethod
    def from_file(cls, ifile, verbose=True, chk_version=True, **kwargs):
        """
        Instantiate the object from an extension in the specified fits file.

        This is a convenience wrapper for :func:`from_hdu`.
        
        Args:
            ifile (:obj:`str`, `Path`_):
                Fits file with the data to read
            verbose (:obj:`bool`, optional):
                Print informational messages
            chk_version (:obj:`bool`, optional):
                Passed to :func:`from_hdu`.
            kwargs (:obj:`dict`, optional):
                Arguments passed directly to :func:`from_hdu`.

        Raises:
            FileNotFoundError:
                Raised if the specified file does not exist.
        """
        _ifile = Path(ifile).resolve()
        if not _ifile.exists():
            raise FileNotFoundError(f'{_ifile} does not exist!')

        if verbose:
            msgs.info(f'Loading {cls.__name__} from {_ifile}')

        # Do it
        with io.fits_open(_ifile) as hdu:
            return cls.from_hdu(hdu, chk_version=chk_version, **kwargs)

    def __repr__(self):
        """ Over-ride print representation

        Returns:
            str: Basics of the Data Container
        """
        repr = '<{:s}: '.format(self.__class__.__name__)
        # Image
        rdict = {}
        for attr in self.datamodel.keys():
            if hasattr(self, attr) and getattr(self, attr) is not None:
                rdict[attr] = True
            else:
                rdict[attr] = False
        repr += ' items={}'.format(rdict)
        repr = repr + '>'
        return repr


def obj_is_data_container(obj):
    """
    Simple method to check whether an object is a data container

    Args:
        obj:

    Returns:
        bool:  True if it is

    """
    return inspect.isclass(obj) and issubclass(obj, DataContainer)



