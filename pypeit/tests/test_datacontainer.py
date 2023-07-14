"""
Module to run tests on datamodel.DataContainer
"""
import sys
import io
import os
import inspect

from IPython import embed

import pytest

import numpy as np

#import pypeit

from astropy.table import Table

from pypeit.datamodel import DataContainer
from pypeit.io import fits_open

#-----------------------------------------------------------------------
# Example derived classes
class BasicContainer(DataContainer):
    version = '1.0.0'
    datamodel = {'vec1': dict(otype=np.ndarray, atype=float, descr='Test'),
                 'meta1': dict(otype=str, decr='test'),
                 'arr1': dict(otype=np.ndarray, atype=float, descr='test')}
    hdu_prefix = 'TST_'

    def __init__(self, vec1, meta1, arr1):
        # All arguments are passed directly to the container
        # instantiation
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        super(BasicContainer, self).__init__({k: values[k] for k in args[1:]}) 

    def _bundle(self):
        # Specify the extension
        return super(BasicContainer, self)._bundle(ext='basic')


class MixedCaseContainer(DataContainer):
    version = '1.0.0'
    datamodel = {'lowercase': dict(otype=np.ndarray, atype=np.integer, descr='Test'),
                 'UPPERCASE': dict(otype=int, decr='test'),
                 'CamelCase': dict(otype=float, decr='test')}

    def __init__(self, lowercase, UPPERCASE, CamelCase):
        # All arguments are passed directly to the container
        # instantiation
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        super(MixedCaseContainer, self).__init__({k: values[k] for k in args[1:]}) 

    def _bundle(self):
        # Specify the extension
        return super(MixedCaseContainer, self)._bundle(ext='mixedcase')


class ImageContainer(DataContainer):
    version = '1.0.0'
    datamodel = {'img1': dict(otype=np.ndarray, atype=float, descr='Test'),
                 'img1_key': dict(otype=str, descr='test'),
                 'img2': dict(otype=np.ndarray, descr='test')}

    def __init__(self, img1, img2, img1_key=None):
        # All arguments are passed directly to the container
        # instantiation
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        super(ImageContainer, self).__init__({k: values[k] for k in args[1:]}) 

    def _bundle(self):
        # img1 and key are put in the first extension, img2 in the
        # second one
        return [ {'img1_key':self['img1_key'], 'img1': self['img1']}, {'img2':self['img2']} ]


class GoodTableContainer(DataContainer):
    version = '1.0.0'
    datamodel = {'tab1': dict(otype=Table, descr='Test'),
                 'tab1len': dict(otype=int, descr='test'),
                 'tab2': dict(otype=Table, descr='test'),
                 'tab2len': dict(otype=int, descr='test')}

    def __init__(self, tab1, tab2):
        # All arguments are passed directly to the container
        # instantiation, but the list is incomplete
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        super(GoodTableContainer, self).__init__({k: values[k] for k in args[1:]}) 

    def _validate(self):
        # Complete the instantiation
        self.tab1len = len(self.tab1)
        self.tab2len = len(self.tab2)

    def _bundle(self):
        # Bundle so there's only one Table per extension
        return [{'tab1len': self.tab1len, 'tab1': self.tab1},
                {'tab2len': self.tab2len, 'tab2': self.tab2}]


class BadTableContainer(GoodTableContainer):
    def _bundle(self):
        # Use default _bundle method, which will try to put both tables
        # in the same extension.  NOTE: Can't use super here because
        # GoodTableContainer doesn't have an 'ext' argument
        return DataContainer._bundle(self, ext='bad')


class GoodMixedTypeContainer(DataContainer):
    version = '1.0.0'
    datamodel = {'tab1': dict(otype=Table, descr='Test'),
                 'tab1len': dict(otype=int, descr='test'),
                 'arr1': dict(otype=np.ndarray, atype=np.integer, descr='test'),
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


class BadInitContainer(DataContainer):
    version = '1.0.0'
    datamodel = {'inp1': dict(otype=np.ndarray, descr='Test'),
                 'inp2': dict(otype=np.ndarray, descr='test'),
                 'out': dict(otype=np.ndarray, descr='test'),
                 'alt': dict(otype=np.ndarray, descr='test')}

    def __init__(self, inp1, inp2, func='add'):
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        super(BadInitContainer, self).__init__({k: values[k] for k in args[1:]}) 


class DubiousInitContainer(DataContainer):
    version = '1.0.0'
    datamodel = {'inp1': dict(otype=np.ndarray, atype=np.integer, descr='Test'),
                 'inp2': dict(otype=np.ndarray, atype=np.integer, descr='test'),
                 'out': dict(otype=np.ndarray, atype=np.integer, descr='test'),
                 'alt': dict(otype=np.ndarray, atype=np.integer, descr='test')}

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

    def _init_internals(self):
        # Because func isn't part of the data model, it won't be part of
        # self if the object is instantiated from a file.  So I have to
        # add it here.  But I don't know what the value of the attribute
        # was for the original object that was written to disk.
        if not hasattr(self, 'func'):
            self.func = None

    def _validate(self):
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
    version = '1.0.0'
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


def test_single_element_array():
    data = BasicContainer(np.arange(1).astype(float), 'length=10', np.arange(10).astype(float))
    hdu = data.to_hdu()
    assert isinstance(hdu[0].data, np.ndarray)
    #
    _data = BasicContainer.from_hdu(hdu[0])
    assert isinstance(_data.vec1, np.ndarray)


def test_basic():

    # Remove file if it exists
    ofile = 'test.fits.gz'
    if os.path.isfile(ofile):
        os.remove(ofile)

    # Instantiate such that number of data table rows would be the same
    # (10)
    data = BasicContainer(np.arange(10).astype(float), 'length=10', np.arange(30).astype(float).reshape(10,3))

    # Instantiation and access tests
    assert list(data.keys()) == ['vec1', 'meta1', 'arr1'], 'Bad keys'
    assert np.array_equal(data.vec1, np.arange(10)), 'Bad instantiation'
    assert np.array_equal(data['vec1'], np.arange(10)), 'Bad instantiation'
    assert np.array_equal(data['arr1'], np.arange(30).reshape(10,3)), 'Bad instantiation'
    assert data.meta1 == 'length=10', 'Bad instantiation'
    assert data.version == BasicContainer.version, 'Bad version'
    # The version number cannot be changed!
    with pytest.raises(TypeError):
        data.version = '2.0'

    # Cannot add elements that aren't part of the datamodel
    with pytest.raises(KeyError):
        data['newvec'] = np.arange(10)
    with pytest.raises(AttributeError):
        data.newvec = np.arange(10)
    # Attempting to access a bad attribute raises an AttributeError
    with pytest.raises(AttributeError):
        test = data.newvec
    # Attempting to access a bad item raises a KeyError
    with pytest.raises(KeyError):
        test = data['newvec']

    # Strict types enforced
    with pytest.raises(TypeError):
        data.vec1 = 3
    with pytest.raises(TypeError):
        data.meta1 = 4.

    # Write to a file
    data.to_file(ofile)

    # Test written data against input
    with fits_open(ofile) as hdu:
        assert len(hdu) == 2, 'Should be two extensions'
        # This tests hdu_prefix
        assert hdu[1].name == 'TST_BASIC', 'Incorrect extension Name'
        assert len(hdu[1].data) == 10, 'Incorrect number of rows'
        assert hdu[1].columns.names == ['vec1', 'arr1'], 'Incorrect column names'
        assert 'meta1' in hdu[1].header, 'Missing header keyword'
        assert hdu[1].header['meta1'] == 'length=10', 'Incorrect header keyword'
   
    # Remove it
    os.remove(ofile)

    # Rows of arrays mismatch, so data is stuffed into a single table
    # row
    data = BasicContainer(np.arange(10).astype(float), 'length=1', np.arange(30).astype(float).reshape(3,10))
    data.to_file(ofile)

    with fits_open(ofile) as hdu:
        assert len(hdu) == 2, 'Should be two extensions'
        assert hdu[1].name == 'TST_BASIC', 'Incorrect extension Name'
        assert len(hdu[1].data) == 1, 'Incorrect number of rows'
        assert hdu[1].columns.names == ['vec1', 'arr1'], 'Incorrect column names'
        assert 'meta1' in hdu[1].header, 'Missing header keyword'
        assert hdu[1].header['meta1'] == 'length=1', 'Incorrect header keyword'

    # Test file input
    _data = BasicContainer.from_file(ofile)
    assert np.array_equal(data.vec1, _data.vec1), 'Bad read'
    assert np.array_equal(data['arr1'], _data['arr1']), 'Bad read'
    assert data.meta1 == _data['meta1'], 'Bad read'
    
    # Clean up
    os.remove(ofile)


def test_case():
    # Remove file if it exists
    ofile = 'test.fits.gz'
    if os.path.isfile(ofile):
        os.remove(ofile)

    # Instantiate
    data = MixedCaseContainer(np.arange(10), 5, 8.2)

    # Attributes and items are case-sensitive
    with pytest.raises(AttributeError):
        data.uppercase = 3
    with pytest.raises(KeyError):
        data['LOWERCASE'] = 3
    with pytest.raises(AttributeError):
        data.cAMELcASE = 3

    # HDU items are not
    data.to_file(ofile)
    with fits_open(ofile) as hdu:
        assert hdu['mixedcase'].header['UpPeRcAsE'] == 5, 'Bad header access'
        assert hdu['mixedcase'].header['CaMeLcAsE'] == 8.2, 'Bad header access'

    # And everything is as it should be when you read it back in
    _data = MixedCaseContainer.from_file(ofile)
    assert _data.UPPERCASE == 5, 'Bad file read'

    # Clean up
    os.remove(ofile)


def test_image():

    img = ImageContainer(np.arange(100).astype(float).reshape(10,10), np.arange(25).reshape(5,5), img1_key='test')
    hdu = img.to_hdu(add_primary=True)
    _img = ImageContainer.from_hdu(hdu)

    assert img.__dict__.keys() == _img.__dict__.keys(), 'Keys are different'
    assert len(hdu) == 3, 'Should be 3 extensions.'

    # No keyword defined
    img = ImageContainer(np.arange(100).astype(float).reshape(10,10), np.arange(25).reshape(5,5))
    assert img.img1_key is None, 'Bad instantiation'
    hdu = img.to_hdu(add_primary=True)
    _img = ImageContainer.from_hdu(hdu)
    assert _img.img1_key is None, 'Bad read'


def test_table():

    x = np.arange(10)
    y = np.arange(10)+5
    z = np.arange(30).reshape(10,3)

    a = np.arange(15).reshape(3,5)
    b = np.zeros(3)
    c = np.full((3,3,3), -1)

    tab1 = Table(data=({'x':x,'y':y,'z':z}), meta={'test':'this'})
    tab2 = Table(data=({'a':a,'b':b,'c':c}), meta={'test':'that'})

    # Check the instantiation
    data = GoodTableContainer(tab1, tab2)
    assert data.tab1.keys() == ['x', 'y', 'z'], 'Bad table save'
    assert data.tab1.meta['test'] == 'this', 'Bad table meta save'

    # Write it to an HDUList
    hdu = data.to_hdu(add_primary=True)
    assert len(hdu) == 3, 'Should be 3 extensions'
    assert [h.name for h in hdu] == ['PRIMARY', 'TAB1', 'TAB2'], 'Bad extension names'
    assert len(hdu['TAB1'].data) == data.tab1len, 'Bad table length'
    assert len(hdu['TAB2'].data) == data.tab2len, 'Bad table length'
    assert hdu['TAB2'].data['c'].shape == (3,3,3), 'Bad column shape'
    assert hdu['TAB2'].header['tab2len'] == 3, 'Bad header data'
    assert hdu['TAB1'].header['TAB1LEN'] == 10, 'Bad header data'

    # Table metadata is also written to the header, which can be
    # accessed with case-insensitive keys
    assert hdu['TAB1'].header['TEST'] == 'this', 'Bad table meta data save'
    assert hdu['TAB1'].header['TesT'] == 'this', 'Bad table meta data save'
    # But the keyword case gets mangled when you read it back in.  This
    # has to do with the Table.read method and I'm not sure there's
    # anything we can do about it without the help of the astropy folks
    _data = GoodTableContainer.from_hdu(hdu)
    with pytest.raises(KeyError):
        _data.tab1.meta['test'] == 'this'
    assert _data.tab1.meta['TEST'] == 'this', 'meta data keywords are re-read in all caps'

    # The BadTableContainer will instantiate fine
    data = BadTableContainer(tab1, tab2)
    assert data.tab1.keys() == ['x', 'y', 'z'], 'Bad table save'
    assert data.tab1.meta['test'] == 'this', 'Bad table meta save'
    # But then it will barf when you try to reformat/write the data
    # because you can't write more than one table to a single HDU
    with pytest.raises(ValueError):
        hdu = data.to_hdu()

def test_mixed():

    x = np.arange(10)
    y = np.arange(10)+5
    z = np.arange(30).reshape(10,3)

    arr1 = np.full((3,3,3), -1)

    tab1 = Table(data=({'x':x,'y':y,'z':z}), meta={'test':'this'})

    # Check the instantiation
    data = GoodMixedTypeContainer(tab1, arr1)
    assert data.tab1.keys() == ['x', 'y', 'z'], 'Bad table save'
    assert data.tab1.meta['test'] == 'this', 'Bad table meta save'
    assert data.arr1shape == (3,3,3), 'Bad shape'

    # Write it to an HDUList
    hdu = data.to_hdu(add_primary=True)
    assert len(hdu) == 3, 'Should be 3 extensions'
    assert [h.name for h in hdu] == ['PRIMARY', 'TAB1', 'ARR1'], 'Bad extension names'
    assert len(hdu['TAB1'].data) == data.tab1len, 'Bad table length'
    assert hdu['ARR1'].shape == data.arr1shape, 'Bad array shape'
    assert hdu['TAB1'].header['TAB1LEN'] == 10, 'Bad header data'
    # Tuples get converted to strings
    assert hdu['ARR1'].header['ARR1SHAPE'] == str((3,3,3)), 'Bad array shape'

    # And they get converted back when they're read in.
    _data = GoodMixedTypeContainer.from_hdu(hdu)
    assert _data.arr1shape == (3,3,3), 'Bad array shape'

    # The BadMixedTypeContainer will instantiate fine
    data = BadMixedTypeContainer(tab1, arr1)
    assert data.tab1.keys() == ['x', 'y', 'z'], 'Bad table save'
    assert data.tab1.meta['test'] == 'this', 'Bad table meta save'
    assert data.arr1shape == (3,3,3), 'Bad shape'
    # But then it will barf when you try to reformat/write the data
    # because you can't write both a Table and an array to a single HDU
    with pytest.raises(ValueError):
        hdu = data.to_hdu()

def test_init():

    x = np.arange(10)
    y = np.arange(10)+5

    # Instantiation of the BadInitContainer should fail because the init
    # arguments all need to be part of the datamodel for it to work.
    with pytest.raises(AttributeError):
        data = BadInitContainer(x,y)

    # This instantiation is fine because DubiousInitContainer handles
    # the fact that some of the arguments to __init__ are not part of
    # the datamodel.
    data = DubiousInitContainer(x,y)
    assert np.array_equal(data.out, data.inp1+data.inp2), 'Bad init'
    # One component of the data model wasn't instantiated, so it will be
    # None
    assert data.alt is None, 'Bad empty datamodel component'

    # Use the other option
    data = DubiousInitContainer(x,y,func='sub')
    assert np.array_equal(data.out, np.full(10,-5)), 'Bad init'

    # Note that this instantiation will set ``out`` because it's None on
    # input.
    stdout = sys.stdout
    sys.stdout = io.StringIO()
    data = DubiousInitContainer(x,y)
    output = sys.stdout.getvalue()
    sys.stdout = stdout
    assert output.strip() == 'Assigning out!', 'Bad instantiation'

    # But if we write it out and read it back in again, the
    # initialization of out will get skipped because it's been read in
    # from the output file.
    hdu = data.to_hdu(add_primary=True)
    stdout = sys.stdout
    sys.stdout = io.StringIO()
    _data = DubiousInitContainer.from_hdu(hdu)
    output = sys.stdout.getvalue()
    sys.stdout = stdout
    assert output.strip() != 'Assigning out!', 'Bad instantiation'

    # The nuance here is that if you instantiate it, change its value
    # away from what it's internals suggest it should be, and write it,
    # it won't be what you expect when you read it back in.  This is a
    # trivial example, but it just points out that while we're strict
    # about the datamodel in these containers, we're not strict about
    # the consistency of the datamodel with respect to the internal
    # operations within the object...
    data = DubiousInitContainer(x,y)
    data.out *= 100
    data.alt = np.arange(20)
    hdu = data.to_hdu(add_primary=True)
    assert not np.array_equal(DubiousInitContainer.from_hdu(hdu).out,
                              DubiousInitContainer(x,y).out), 'Bad read'

    # The reason the DubiousInitContainer is called that way is because
    # it has attributes that cannot be reinstantiated from what's
    # written to disk, meaning these DataContainers aren't really fully
    # formed objects unless all of its relevant attributes are
    # components of the data model.
    _data = DubiousInitContainer.from_hdu(hdu)
    with pytest.raises(AssertionError):
        assert _data.func == DubiousInitContainer(x,y).func

    # This is solved by adding func to the datamodel
    data = ComplexInitContainer(x,y)
    _data = ComplexInitContainer.from_hdu(data.to_hdu(add_primary=True))
    assert data.func == _data.func, 'Bad read'

