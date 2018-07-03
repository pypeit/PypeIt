# encoding: utf-8
"""
Defines parameter sets used to set the behavior for core pypit
functionality.

There are two main parameter set groups, those used for reducing the
data:

    - RunPar: Parameters specific to a given execution of PypIt

    - ReducePar: Parameters the define how PypIt should perform the
      reductions.  These include general parameters and the following
      parameter subsets:
        - OverscanPar: Methods used for the overscan subtraction.
        - FlatFieldPar: Methods used for field-flattening.
        - SkySubtractionPar: Methods used for sky subtraction.
        - FlexurePar: Methods used for flexure correction.
        - WavelengthCalibrationPar: Methods used for constructing the
          wavelength solution.
        - FluxCalibrationPar: Methods used for flux calibration.

    - FrameGroupPar: Sets parameters that are used to group and combined
      frames.
        - CombineFramesPar - Parameters for combining frames

    - WavelengthSolutionPar: Parameters used for constructing the wavelength
      solution

    - TraceSlitsPar: Parameters for tracing slits
        - PCAPar: PCA parameters

    - TraceTiltsPar: Parameters for tracing tilts of arc lines

    - TraceObjectsPar: Parameters for tracing objects

    - ExtractObjectsPar: Parameters for extracting 1D object spectra
        - ManualExtractionPar: Parameters needed for manual extraction
          of 1D spectra

And those that define the parameters of the instrument used to collect
the data:

    - DetectorPar: Specifies the properties of a given detector

    - InstrumentPar: Specifies the instrument used for a given
      observation and contains the list of detectors for the camera.

    - FrameFitsPar: Defines the set of fits files properties that all
      fits files from the given instrument should have.

    - FrameIDPar: Defines a set of parameters used to identify fits
      files as a certain frame type.

These are collected into the main, high-level paramter set, called
PypitPar.


**New Parameters**:

To add a new parameter, let's call it `foo`, to any of the provided
parameter sets:
    - Add `foo=None` to the __init__ method of the relevant parameter
      set.

    - Add any default value (the default value is None unless you set
      it), options list, data type, and description to the body of the
      __init__ method.  Like so::

        defaults['foo'] = 'bar'
        options['foo'] = [ 'bar', 'boo' ]
        dtypes['foo'] = str
        descr['foo'] = 'foo? who you callin a foo!  ' \
                       'Options are: {0}'.format(', '.join(options['foo']))

    - Add the parameter to the `from_dict` method:
    
        - If the parameter is something that does not require
          instantiation, add the keyword to the `parkeys` list in the
          `from_dict` method

        - If the parameter is another ParSet or requires instantiation,
          provide the instantiation.  For example, see how the
          `CombineFramesPar` parameter set is defined in the
          `FrameGroupPar` class.  E.g.::

            pk = 'foo'
            kwargs[pk] = FooPar.from_dict(cfg[pk]) if pk in k else None

**New Parameter Sets:**

To add an entirely new parameter set, use one of the existing parameter
sets as a template, then add the parameter set to `PypitPar` (assuming
you want it to be accessed throughout the code.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

try:
    basestring
except NameError:
    basestring = str

import os
import glob
import warnings
from pkg_resources import resource_filename
import inspect

from collections import OrderedDict

import numpy

from configobj import ConfigObj
from astropy.time import Time

from .parset import ParSet

#-----------------------------------------------------------------------------
# Helper functions

def _pypit_root_directory():
    """
    Get the root directory for the PYPIT source distribution.

    .. todo::
        - Set this in __init__.py

    Returns:
        str: Root directory to PYPIT

    Raises:
        OSError: Raised if `pkg_resources.resource_filename` fails.
    """
    try:
        # Get the directory with the pypit source code
        code_dir = resource_filename('pypit', '')
    except:
        # TODO: pypit should always be installed as a package, so is
        # this try/except block necessary?
        raise OSError('Could not find PYPIT package!')
    # Root directory is one level up from source code
    return os.path.split(code_dir)[0]
#    for p in sys.path:
#        if 'PYPIT'.lower() in p.lower() \
#                and len(glob.glob(os.path.join(p, 'pypit', 'pypit.py'))) == 1:
#            return p
#    raise OSError('Could not find PYPIT in system path.')


def _eval_ignore():
    """Provides a list of strings that should not be evaluated."""
    return [ 'open', 'file', 'dict' ]


def _recursive_dict_evaluate(d):
    """
    Recursively run :func:`eval` on each element of the provided
    dictionary.

    A raw read of a configuration file with `ConfigObj` results in a
    dictionary that contains strings or lists of strings.  However, when
    assigning the values for the various ParSets, the `from_dict`
    methods expect the dictionary values to have the appropriate type.
    E.g., the ConfigObj will have something like d['foo'] = '1', when
    the `from_dict` method expects the value to be an integer (d['foo']
    = 1).

    This function tries to evaluate *all* dictionary values, except for
    those listed above in the :func:`_eval_ignore` function.  Any value
    in this list or where::

        eval(d[k]) for k in d.keys()

    raises an exception is returned as the original string.

    This is currently only used in :func:`PypitPar.from_cfg_file`; see
    further comments there.

    Args:
        d (dict):
            Dictionary of values to evaluate

    Returns:
        dict: Identical to input dictionary, but with all string values
        replaced with the result of `eval(d[k])` for all `k` in
        `d.keys()`.
    """
    ignore = _eval_ignore()
    for k in d.keys():
        if isinstance(d[k], dict):
           d[k] = _recursive_dict_evaluate(d[k])
        elif isinstance(d[k], list):
            replacement = []
            for v in d[k]:
                if v in ignore:
                    replacement += [ v ]
                else:
                    try:
                        replacement += [ eval(v) ]
                    except:
                        replacement += [ v ]
            d[k] = replacement
        else:
            try:
                d[k] = eval(d[k]) if d[k] not in ignore else d[k]
            except:
                pass

    return d


#def _recursive_dict_unicode2str(d):
#    """
#    Recursively convert any unicode to a string.
#    """
#    for k in d.keys():
#        if isinstance(d[k], dict):
#           d[k] = _recursive_dict_unicode2str(d[k])
#        elif isinstance(d[k], list):
#            replacement = []
#            for v in d[k]:
#                if isinstance(v, unicode):
#                    replacement += [ str(v) ]
#                    continue
#                replacement += [ v ]
#            d[k] = replacement
#        else:
#            if isinstance(d[k], unicode):
#                d[k] = str(d[k])
#    return d


def _get_parset_list(cfg, pk, parsetclass):
    """
    Create a list of ParSets based on a root keyword for a set of
    defined groups in the configuration file.
    
    For example, the :class:`InstrumentPar` group allows for a list of
    detectors (:class:`DetectorPar`) with keywords like `detector1`,
    `detector2`, etc.  This function parses the provided configuration
    object (`cfg`) to find any sections with `detector` (`pk`) as its
    root.  The remainder of the section name must be able to be
    converted to an integer and the section itself must be able to setup
    an instance of `parsetclass`.  The sections must be number
    sequentially from 1..N.  E.g., the :class:`InstrumentPar`
    configuration file cannot have `dectector1` and `detector3`, but no
    `detector2`.  The call to setup the detectors in the
    :class:`InstrumentPar` is::

        kwargs['detector'] = _get_parset_list(cfg, 'detector', DetectorPar)

    Args:
        cfg (:class:`ConfigObj`, :obj:`dict`):
            The top-level configuration that defines a list of
            sub-ParSets.
        pk (str):
            The root of the keywords used to set a list of sub-ParSets.
        parsetclass (:class:`pypit.par.parset.ParSet`):
            The class used to construct each element in the list of
            parameter subsets.  The class **must** have a `from_dict`
            method that instantiates the
            :class:`pypit.par.parset.ParSet` based on the provide
            subsection/subdict from cfg.

    Returns:
        list: A list of instances of `parsetclass` parsed from the
        provided configuration data.

    Raises:
        ValueError:
            Raised if the indices of the subsections are not sequential
            and 1-indexed.
    """
    # Get the full list of keys
    k = cfg.keys()

    # Iterate through the list of keys to find the appropriate sub
    # parameter sets and their order.
    par = []
    order = []
    for _k in k:
        if _k == pk and cfg[_k] is None:
            continue
        if pk in _k:
            try:
                # Get the order for this subgroup (e.g., 2 for
                # 'detector2'
                order += [ int(_k.replace(pk,'')) ]
                # And instantiate the parameter set
                par += [ parsetclass.from_dict(cfg[_k]) ]
            except:
                continue

    if len(par) > 0:
        # Make sure the instances are correctly sorted and sequential
        srt = numpy.argsort(order)
        if numpy.any(numpy.array(order)[srt]-1 != numpy.arange(order[srt[-1]])):
            raise ValueError('Parameter set series must be sequential and 1-indexed.')
        # Return the sorted instances
        return [par[i] for i in srt]

    # No such subsets were defined, so return a null result
    return None

#-----------------------------------------------------------------------------
# Reduction ParSets

class FrameGroupPar(ParSet):
    def __init__(self, frametype=None, useframe=None, number=None, combine=None):
        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['frametype'] = 'bias'
        options['frametype'] = FrameGroupPar.valid_frame_types()
        dtypes['frametype'] = basestring
        descr['frametype'] = 'Frame type.  ' \
                             'Options are: {0}'.format(', '.join(options['frametype']))

        dtypes['useframe'] = basestring
        descr['useframe'] = 'A master calibrations file to use if it exists'

        defaults['number'] = 0
        dtypes['number'] = int
        descr['number'] = 'Number of frames to use of this type'

        defaults['combine'] = CombineFramesPar()
        dtypes['combine'] = [ParSet, dict]
        descr['combine'] = 'Parameters used when combining frames of this type'

        # Instantiate the parameter set
        super(FrameGroupPar, self).__init__(list(pars.keys()),
                                            values=list(pars.values()),
                                            defaults=list(defaults.values()),
                                            options=list(options.values()),
                                            dtypes=list(dtypes.values()),
                                            descr=list(descr.values()))

        self.validate()

    @classmethod
    def from_dict(cls, frametype, cfg):
        k = cfg.keys()
        parkeys = [ 'useframe', 'number' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        pk = 'combine'
        kwargs[pk] = CombineFramesPar.from_dict(cfg[pk]) if pk in k else None
        return cls(frametype=frametype, **kwargs)

    @staticmethod
    def valid_frame_types():
        """
        Return the list of valid frame types.
        """
        return [ 'bias', 'pixelflat', 'arc', 'pinhole', 'trace', 'standard', 'science' ]

    def validate(self):
        pass


class CombineFramesPar(ParSet):
    def __init__(self, match=None, method=None, satpix=None, cosmics=None, n_lohi=None,
                 sig_lohi=None, replace=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['match'] = -1
        dtypes['match'] = [int, float]
        descr['match'] = 'Match frames with pixel counts that are within N-sigma of one ' \
                         'another, where match=N below.  If N < 0, nothing is matched.'

        defaults['method'] = 'mean'
        options['method'] = CombineFramesPar.valid_methods()
        dtypes['method'] = basestring
        descr['method'] = 'Method used to combine frames.  Options are: {0}'.format(
                                       ', '.join(options['method']))

        defaults['satpix'] = 'reject'
        options['satpix'] = CombineFramesPar.valid_saturation_handling()
        dtypes['satpix'] = basestring
        descr['satpix'] = 'Handling of saturated pixels.  Options are: {0}'.format(
                                       ', '.join(options['satpix']))

        defaults['cosmics'] = 20.0
        dtypes['cosmics'] = [int, float]
        descr['cosmics'] = 'Sigma level to reject cosmic rays (<= 0.0 means no CR removal)'

        # TODO: Test the list length...
        defaults['n_lohi'] = [0, 0]
        dtypes['n_lohi'] = list
        descr['n_lohi'] = 'Number of pixels to reject at the lowest and highest ends of the ' \
                          'distribution; i.e., n_lohi = low, high.  Use None for no limit.'

        defaults['sig_lohi'] = [3.0, 3.0]
        dtypes['sig_lohi'] = list
        descr['sig_lohi'] = 'Sigma-clipping level at the low and high ends of the distribution; ' \
                            'i.e., sig_lohi = low, high.  Use None for no limit.'

        defaults['replace'] = 'maxnonsat'
        options['replace'] = CombineFramesPar.valid_rejection_replacements()
        dtypes['replace'] = basestring
        descr['replace'] = 'If all pixels are rejected, replace them using this method.  ' \
                           'Options are: {0}'.format(', '.join(options['replace']))

        # Instantiate the parameter set
        super(CombineFramesPar, self).__init__(list(pars.keys()),
                                               values=list(pars.values()),
                                               defaults=list(defaults.values()),
                                               options=list(options.values()),
                                               dtypes=list(dtypes.values()),
                                               descr=list(descr.values()))
        
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'match', 'method', 'satpix', 'cosmics', 'n_lohi', 'sig_lohi', 'replace' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_methods():
        """
        Return the valid methods for combining frames.
        """
        return [ 'mean', 'median', 'weightmean' ]

    @staticmethod
    def valid_saturation_handling():
        """
        Return the valid approachs to handling saturated pixels.
        """
        return [ 'reject', 'force', 'nothing' ]

    @staticmethod
    def valid_rejection_replacements():
        """
        Return the valid replacement methods for rejected pixels.
        """
        return [ 'min', 'max', 'mean', 'median', 'weightmean', 'maxnonsat' ]

    def validate(self):
        pass


class OverscanPar(ParSet):
    def __init__(self, method=None, params=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['method'] = 'savgol'
        options['method'] = OverscanPar.valid_methods()
        dtypes['method'] = basestring
        descr['method'] = 'Method used to fit the overscan.  ' \
                          'Options are: {0}'.format(', '.join(options['method']))
        
        defaults['params'] = [5, 65]
        dtypes['params'] = [int, list]
        descr['params'] = 'Parameters for the overscan subtraction.  For \'polynomial\', set ' \
                          'params = order, number of pixels, number of repeats ; for ' \
                          '\'savgol\', set params = order, window size ; for \'median\', set ' \
                          'params = None or omit the keyword.'

        # Instantiate the parameter set
        super(OverscanPar, self).__init__(list(pars.keys()),
                                          values=list(pars.values()),
                                          defaults=list(defaults.values()),
                                          options=list(options.values()),
                                          dtypes=list(dtypes.values()),
                                          descr=list(descr.values()))

        # Check the parameters match the method requirements
        self.validate()


    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'method', 'params' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_methods():
        """
        Return the valid overscane methods.
        """
        return [ 'polynomial', 'savgol', 'median' ]

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['method'] is None:
            return
        if self.data['params'] is None:
            raise ValueError('No overscan method parameters defined!')

        if self.data['method'] == 'polynomial' and len(self.data['params']) != 3:
            raise ValueError('For polynomial overscan method, set params = order, number of '
                             'pixels, number of repeats')

        if self.data['method'] == 'savgol' and len(self.data['params']) != 2:
            raise ValueError('For savgol overscan method, set params = order, window size')
            
        if self.data['method'] == 'median' and self.data['params'] is not None:
            warnings.warn('No parameters necessary for median overscan method.  Ignoring input.')


class FlatFieldPar(ParSet):
    def __init__(self, frame=None, method=None, params=None, twodpca=None):
    
        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        # TODO: Provide a list of valid masters to use as options?
        defaults['frame'] = 'pixelflat'
        dtypes['frame'] = basestring
        descr['frame'] = 'Frame to use for field flattening.  Options are: pixelflat, pinhole, ' \
                         'or a specified master calibration file.'

        defaults['method'] = 'bspline'
        options['method'] = FlatFieldPar.valid_methods()
        dtypes['method'] = basestring
        descr['method'] = 'Method used to flat field the data; use None to skip flat-fielding.  ' \
                          'Options are: None, {0}'.format(', '.join(options['method']))

        defaults['params'] = 20
        dtypes['params'] = [int, list]
        descr['params'] = 'Flat-field method parameters.  For \'PolyScan\', set params = order, ' \
                          'numPixels, repeat ; for bspline, set params = spacing '

        # TODO:  How is twodpca used?  Is it just another method that is
        # only used for ARMED?  Could we remove twodpca and just add
        # another method option?
        defaults['twodpca'] = 0
        dtypes['twodpca'] = int
        descr['twodpca'] = 'Perform a simple 2D PCA on the echelle blaze fits if the value of ' \
                           'this argument is >1. The argument value is equal to the number of ' \
                           'PCA components. 0 means that no PCA will be performed.  **This is ' \
                           'only used with ARMED pipeline.'

        # Instantiate the parameter set
        super(FlatFieldPar, self).__init__(list(pars.keys()),
                                           values=list(pars.values()),
                                           defaults=list(defaults.values()),
                                           options=list(options.values()),
                                           dtypes=list(dtypes.values()),
                                           descr=list(descr.values()))

        # Check the parameters match the method requirements
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'frame', 'method', 'params', 'twodpca' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_frames():
        """
        Return the valid frame types.
        """
        return ['pixelflat', 'pinhole']

    @staticmethod
    def valid_methods():
        """
        Return the valid flat-field methods
        """
        return [ 'PolyScan', 'bspline' ]

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['method'] == 'PolyScan' and len(self.data['params']) != 3:
            raise ValueError('For PolyScan method, set params = order, number of '
                             'pixels, number of repeats')

        if self.data['method'] == 'bspline' and not isinstance(self.data['params'], int):
            raise ValueError('For bspline method, set params = spacing (integer).')
            
        if self.data['frame'] in FlatFieldPar.valid_frames() or self.data['frame'] is None:
            return

        if not os.path.isfile(self.data['frame']):
            raise ValueError('Provided frame file name does not exist: {0}'.format(
                                self.data['frame']))


class FlexurePar(ParSet):
    def __init__(self, method=None, maxshift=None, spectrum=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        defaults['method'] = 'boxcar'
        options['method'] = FlexurePar.valid_methods()
        dtypes['method'] = basestring
        descr['method'] = 'Method used to correct for flexure. Use None for no correction.  If ' \
                          'slitcen is used, the flexure correction is performed before the ' \
                          'extraction of objects.  ' \
                          'Options are: None, {0}'.format(', '.join(options['method']))

        defaults['maxshift'] = 20
        dtypes['maxshift'] = [int, float]
        descr['maxshift'] = 'Maximum allowed flexure shift in pixels.'

        dtypes['spectrum'] = basestring
        descr['spectrum'] = 'Archive sky spectrum to be used for the flexure correction.'

        # Instantiate the parameter set
        super(FlexurePar, self).__init__(list(pars.keys()),
                                         values=list(pars.values()),
                                         defaults=list(defaults.values()),
                                         options=list(options.values()),
                                         dtypes=list(dtypes.values()),
                                         descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'method', 'maxshift', 'spectrum' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_frames():
        """
        Return the valid frame types.
        """
        return ['pixelflat', 'pinhole']

    @staticmethod
    def valid_methods():
        """
        Return the valid flat-field methods
        """
        return [ 'boxcar', 'slitcen' ]

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['spectrum'] is not None and not os.path.isfile(self.data['spectrum']):
            raise ValueError('Provided archive spectrum does not exist: {0}.'.format(
                             self.data['spectrum']))


class WavelengthCalibrationPar(ParSet):
    def __init__(self, medium=None, refframe=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        defaults['medium'] = 'vacuum'
        options['medium'] = WavelengthCalibrationPar.valid_media()
        dtypes['medium'] = basestring
        descr['medium'] = 'Medium used when wavelength calibrating the data.  ' \
                          'Options are: {0}'.format(', '.join(options['medium']))

#        defaults['refframe'] = 'heliocentric'
        options['refframe'] = WavelengthCalibrationPar.valid_reference_frames()
        dtypes['refframe'] = basestring
        descr['refframe'] = 'Frame of reference for the wavelength calibration.  ' \
                            'Options are: {0}'.format(', '.join(options['refframe']))

        # Instantiate the parameter set
        super(WavelengthCalibrationPar, self).__init__(list(pars.keys()),
                                                    values=list(pars.values()),
                                                    defaults=list(defaults.values()),
                                                    options=list(options.values()),
                                                    dtypes=list(dtypes.values()),
                                                    descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'medium', 'refframe' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_media():
        """
        Return the valid flat-field methods
        """
        return [ 'vacuum', 'air' ]

    @staticmethod
    def valid_reference_frames():
        """
        Return the valid frame types.
        """
        return [ 'heliocentric', 'barycentric' ]
    
    def validate(self):
        pass


class FluxCalibrationPar(ParSet):
    def __init__(self, flux=None, nonlinear=None, sensfunc=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['flux'] = False
        dtypes['flux'] = bool
        descr['flux'] = 'Flag to perform flux calibration'

        # TODO: I don't think this is used anywhere
        defaults['nonlinear'] = False
        dtypes['nonlinear'] = bool
        descr['nonlinear'] = 'Perform a non-linear correction.  Requires a series of ' \
                             'pixelflats of the same lamp and setup and with a variety of ' \
                             'exposure times and count rates in every pixel.'

        dtypes['sensfunc'] = basestring
        descr['sensfunc'] = 'YAML file with an existing calibration function'

        # Instantiate the parameter set
        super(FluxCalibrationPar, self).__init__(list(pars.keys()),
                                                 values=list(pars.values()),
                                                 defaults=list(defaults.values()),
                                                 dtypes=list(dtypes.values()),
                                                 descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'flux', 'nonlinear', 'sensfunc' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['sensfunc'] is not None and not os.path.isfile(self.data['sensfunc']):
            raise ValueError('Provided sensitivity function does not exist: {0}.'.format(
                             self.data['sensfunc']))


class SkySubtractionPar(ParSet):
    def __init__(self, method=None, params=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['method'] = 'bspline'
        options['method'] = SkySubtractionPar.valid_methods()
        dtypes['method'] = basestring
        descr['method'] = 'Method used to for sky subtraction.  ' \
                          'Options are: None, {0}'.format(', '.join(options['method']))

        defaults['params'] = 20
        dtypes['params'] = int
        descr['params'] = 'Sky-subtraction method parameters.  For bspline, set params = spacing.'

        # Instantiate the parameter set
        super(SkySubtractionPar, self).__init__(list(pars.keys()),
                                                values=list(pars.values()),
                                                defaults=list(defaults.values()),
                                                options=list(options.values()),
                                                dtypes=list(dtypes.values()),
                                                descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'method', 'params' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_methods():
        """
        Return the valid sky-subtraction methods
        """
        return [ 'bspline' ]

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['method'] == 'bspline' and not isinstance(self.data['params'], int):
            raise ValueError('For bspline sky-subtraction method, set params = spacing (integer).')


class PCAPar(ParSet):
    def __init__(self, pcatype=None, params=None, extrapolate=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        # TODO: Change pcatype to spacing='irregular' or
        # spacing='smooth'?
        defaults['pcatype'] = 'pixel'
        options['pcatype'] = PCAPar.valid_types()
        dtypes['pcatype'] = basestring
        descr['pcatype'] = 'Select to perform the PCA using the pixel position (pcatype=pixel) ' \
                           'or by spectral order (pcatype=order).  Pixel positions can be used ' \
                           'for multi-object spectroscopy where the gap between slits is ' \
                           'irregular.  Order is used for echelle spectroscopy or for slits ' \
                           'with separations that are a smooth function of the slit number.'

        defaults['params'] = [ 3, 2, 1, 0, 0, 0 ]
        dtypes['params'] = list
        descr['params'] = 'Order of the polynomials to be used to fit the principle ' \
                          'components.  TODO: Provide more explanation'

        defaults['extrapolate'] = [0, 0]
        dtypes['extrapolate'] = list
        descr['extrapolate'] = 'The number of extra orders to predict in the negative (first ' \
                               'number) and positive (second number) direction.  Must be two ' \
                               'numbers in the list and they must be integers.'

        # Instantiate the parameter set
        super(PCAPar, self).__init__(list(pars.keys()),
                                     values=list(pars.values()),
                                     defaults=list(defaults.values()),
                                     options=list(options.values()),
                                     dtypes=list(dtypes.values()),
                                     descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'pcatype', 'params', 'extrapolate' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_types():
        """
        Return the valid PCA types.
        """
        return ['pixel', 'order']

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        if len(self.data['extrapolate']) != 2:
            raise ValueError('Extrapolate must be a list with two values.')
        for e in self.data['extrapolate']:
            if not isinstance(e, int):
                raise ValueError('Extrapolate values must be integers.')


class ManualExtractionPar(ParSet):
    def __init__(self, frame=None, params=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        dtypes['frame'] = basestring
        descr['frame'] = 'The name of the fits file for a manual extraction'

        dtypes['params'] = list
        descr['params'] = 'Parameters of the manual extraction.  For example, params = ' \
                          '1,1000,500,10,10 specifies the following behavior: 1 is the ' \
                          'detector number, 1000 is the spatial location that the trace must ' \
                          'go through, 500 is the spectral location that the trace must go ' \
                          'through, and the last two numbers (10,10) are the widths around the ' \
                          'stated (spatial,spectral) location that should also be in the trace.'

        # Instantiate the parameter set
        super(ManualExtractionPar, self).__init__(list(pars.keys()),
                                                  values=list(pars.values()),
                                                  dtypes=list(dtypes.values()),
                                                  descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'frame', 'params' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['params'] is not None and len(self.data['params']) != 5:
            raise ValueError('There must be 5 manual extraction parameters.')
        if self.data['frame'] is not None and not os.path.isfile(self.data['frame']):
            raise FileNotFoundError('Manual extraction frame does not exist: {0}'.format(
                                    self.data['frame']))


class RunPar(ParSet):
    """
    Parameters specific to a given execution of PypIt.
    """
    def __init__(self, ncpus=None, calcheck=None, calwin=None, setup=None, qa=None, preponly=None,
                 stopcheck=None, useIDname=None, verbosity=None, caldir=None, scidir=None,
                 qadir=None, sortdir=None, overwrite=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['ncpus'] = 1
        dtypes['ncpus']   = int
        descr['ncpus']    = 'Number of CPUs to use (-1 means all bar one CPU, -2 means all bar ' \
                            'two CPUs)'

        defaults['calcheck'] = False
        dtypes['calcheck']   = bool
        descr['calcheck']    = 'Flag to skip the data reduction and only check that all ' \
                               'calibration data are present'

        defaults['calwin'] = 0
        dtypes['calwin']   = [int, float]
        descr['calwin'] = 'The window of time in hours to search for calibration frames for a ' \
                          'science frame'

        defaults['setup'] = False
        dtypes['setup']   = bool
        descr['setup']    = 'If True, run in setup mode.  Useful to parse files when starting ' \
                            'reduction on a large set of data'

        defaults['qa'] = False
        dtypes['qa']   = bool
        descr['qa']    = 'If True, run quality control in real time.  Otherwise, checks are ' \
                         'saved to disk for later inspection'
        
        defaults['preponly'] = False
        dtypes['preponly']   = bool
        descr['preponly']    = 'If True, pypit will prepare the calibration frames; if false, ' \
                               'will only reduce the science frames'

        defaults['stopcheck'] = False
        dtypes['stopcheck'] = bool
        descr['stopcheck'] = 'If True, pypit will stop and require a user carriage return at ' \
                             'every quality control check'

        defaults['useIDname'] = False
        dtypes['useIDname']   = bool
        descr['useIDname']    = 'If True, file sorting will ensure that the idname is made'

       
        # TODO: Is this needed?  Should be used when instantiating msgs
        defaults['verbosity'] = 2
        dtypes['verbosity']   = int
        descr['verbosity']    = 'Level of screen output: 0 supresses all output; 1 provides ' \
                                'high-level output; 2 provides all output'

        defaults['caldir'] = 'MF'
        dtypes['caldir'] = basestring
        descr['caldir'] = 'Directory relative to calling directory to write master files.'
        
        defaults['scidir'] = 'Science'
        dtypes['scidir'] = basestring
        descr['scidir'] = 'Directory relative to calling directory to write science files.'
        
        defaults['qadir'] = 'QA'
        dtypes['qadir'] = basestring
        descr['qadir'] = 'Directory relative to calling directory to write qa files.'

        defaults['sortdir'] = None
        dtypes['sortdir'] = basestring
        descr['sortdir'] = 'File for the details of the sorted files.  If None, no output is ' \
                           'created'

        defaults['overwrite'] = False
        dtypes['overwrite'] = bool
        descr['overwrite'] = 'Flag to overwrite any existing output files'

        # Instantiate the parameter set
        super(RunPar, self).__init__(list(pars.keys()),
                                     values=list(pars.values()),
                                     defaults=list(defaults.values()),
                                     options=list(options.values()),
                                     dtypes=list(dtypes.values()),
                                     descr=list(descr.values()))
#                                     , cfg_section='run', cfg_comment='Execution options')
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'ncpus', 'calcheck', 'calwin', 'setup', 'qa', 'preponly', 'stopcheck',
                    'useIDname', 'verbosity', 'caldir', 'scidir', 'qadir', 'sortdir', 'overwrite' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        pass


class ReducePar(ParSet):
    """
    Parameters specific to the reduction procedures used by PypIt.
    """
    def __init__(self, spectrograph=None, pipeline=None, detnum=None, masters=None, setup=None,
                 trim=None, badpix=None, slit_center_frame=None, slit_edge_frame=None,
                 overscan=None, flatfield=None, flexure=None, wavecalib=None, fluxcalib=None,
                 skysubtract=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['spectrograph'] = 'KECK_LRISb'
        options['spectrograph'] = ReducePar.valid_spectrographs()
        dtypes['spectrograph'] = basestring
        descr['spectrograph'] = 'Spectrograph that provided the data to be reduced.  ' \
                                'Options are: {0}'.format(', '.join(options['spectrograph']))

        options['pipeline'] = ReducePar.valid_pipelines()
        dtypes['pipeline'] = basestring
        descr['pipeline'] = 'Pipeline options that pypit can use for reductions.  ' \
                            'Options are: {0}'.format(', '.join(options['pipeline']))

        dtypes['detnum'] = int
        descr['detnum'] = 'Restrict reduction to a single detector with this index'

        options['masters'] = ReducePar.allowed_master_options()
        dtypes['masters'] = basestring
        descr['masters'] = 'Treatment of master frames.  Use None to select the default ' \
                           'behavior (which is?), \'reuse\' to use any existing masters, and ' \
                           '\'force\' to __only__ use master frames.  ' \
                           'Options are: None, {0}'.format(', '.join(options['masters']))

        dtypes['setup'] = basestring
        descr['setup'] = 'If masters=\'force\', this is the setup name to be used: e.g., ' \
                         'C_02_aa .  The detector number is ignored but the other information ' \
                         'must match the Master Frames in the master frame folder.'

        defaults['trim'] = True
        dtypes['trim'] = bool
        descr['trim'] = 'Trim the frame to isolate the data'

        defaults['badpix'] = True
        dtypes['badpix'] = bool
        descr['badpix'] = 'Make a bad pixel mask? Bias frames must be provided.'

        # TODO: What are the allowed center and edge trace frames?
        defaults['slit_center_frame'] = 'trace'
        dtypes['slit_center_frame'] = basestring
        descr['slit_center_frame'] = 'The frame that should be used to trace the slit ' \
                                     'centroid.  A master calibrations file can also be ' \
                                     'specified.'

        defaults['slit_edge_frame'] = 'trace'
        dtypes['slit_edge_frame'] = basestring
        descr['slit_edge_frame'] = 'The frame that should be used to trace the slit edges.  ' \
                                   'A master calibrations file can also be specified.'

        defaults['overscan'] = OverscanPar()
        dtypes['overscan'] = [ ParSet, dict ]
        descr['overscan'] = 'Parameters used to fit the overscan region.'

        defaults['flatfield'] = FlatFieldPar()
        dtypes['flatfield'] = [ ParSet, dict ]
        descr['flatfield'] = 'Parameters used to set the flat-field procedure'

        defaults['flexure'] = FlexurePar()
        dtypes['flexure'] = [ ParSet, dict ]
        descr['flexure'] = 'Parameters used to set the flexure-correction procedure'

        defaults['wavecalib'] = WavelengthCalibrationPar()
        dtypes['wavecalib'] = [ ParSet, dict ]
        descr['wavecalib'] = 'Parameters used to set the wavelength calibration to provide for ' \
                             'the output spectra'

        defaults['fluxcalib'] = FluxCalibrationPar()
        dtypes['fluxcalib'] = [ ParSet, dict ]
        descr['fluxcalib'] = 'Parameters used to set the flux-calibration procedure'
        
        defaults['skysubtract'] = SkySubtractionPar()
        dtypes['skysubtract'] = [ ParSet, dict ]
        descr['skysubtract'] = 'Parameters used to set the sky-subtraction procedure'

        # Instantiate the parameter set
        super(ReducePar, self).__init__(list(pars.keys()),
                                        values=list(pars.values()),
                                        defaults=list(defaults.values()),
                                        options=list(options.values()),
                                        dtypes=list(dtypes.values()),
                                        descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()

        # Basic keywords
        parkeys = [ 'spectrograph', 'pipeline', 'detnum', 'masters', 'setup', 'trim', 'badpix',
                    'slit_center_frame', 'slit_edge_frame' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None

        # Keywords that are ParSets
        pk = 'overscan'
        kwargs[pk] = OverscanPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'flatfield'
        kwargs[pk] = FlatFieldPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'flexure'
        kwargs[pk] = FlexurePar.from_dict(cfg[pk]) if pk in k else None
        pk = 'wavecalib'
        kwargs[pk] = WavelengthCalibrationPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'fluxcalib'
        kwargs[pk] = FluxCalibrationPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'skysubtract'
        kwargs[pk] = SkySubtractionPar.from_dict(cfg[pk]) if pk in k else None

        return cls(**kwargs)

    @staticmethod
    def valid_spectrographs():
        """
        Return the list of allowed spectrographs for pypit reductions.
        The valid spectrographs are determined by finding the *.cfg
        files in PYPIT source directory structure, and any *.cfg files
        in the current working directory.
        
        .. todo::
            - Remove default_spectrograph.cfg
        """
        # Find the spectrograph files included in the distribution
        pypit_root = _pypit_root_directory()
        spec_dir = os.path.join(pypit_root, 'pypit', 'config', 'spectrographs')
        cfg_files = glob.glob(os.path.join(spec_dir, '*_spectrograph.cfg'))

        # Find any user spectrograph files
        user_cfg_files = glob.glob(os.path.join('*_spectrograph.cfg'))
        if len(user_cfg_files) > 0:
            warnings.warn('Found *_spectrograph.cfg files in current working directory.  '
                          'Spectrograph options will include these files.')
            cfg_files += user_cfg_files
        
        return [ '_'.join((f.split('/')[-1]).split('_')[:-1]) for f in cfg_files ]

    @staticmethod
    def spectrograph_config_file(key, verbose=False):
        """
        Return the list of allowed spectrographs for pypit reductions.
        The valid spectrographs are determined by finding the *.cfg
        files in PYPIT source directory structure, and any *.cfg files
        in the current working directory.
        
        .. todo::
            - Remove default_spectrograph.cfg
        """
        # Name of the file to find
        file_name = '{0}_spectrograph.cfg'.format(key)

        # First try to find the file in the local directory
        cfg_file = os.path.join(os.getcwd(), file_name)
        if os.path.isfile(cfg_file):
            if verbose:
                print('Found: {0}'.format(cfg_file))
            return cfg_file

        # Then try to find the file in the distribution
        pypit_root = _pypit_root_directory()
        spec_dir = os.path.join(pypit_root, 'pypit', 'config', 'spectrographs')
        cfg_file = os.path.join(spec_dir, file_name)
        if os.path.isfile(cfg_file):
            if verbose:
                print('Found: {0}'.format(cfg_file))
            return cfg_file

        # Could not find the file!
        raise ValueError('Could not find associated configuration file in either the local'
                         'directory or the pypit source distribution for '
                         'spectrograph {0}!'.format(key))

    @staticmethod
    def valid_pipelines():
        """Return the list of allowed pipelines within pypit."""
        return [ 'ARMS', 'ARMED' ]

    @staticmethod
    def allowed_master_options():
        """Return the allowed handling methods for the master frames."""
        return [ 'reuse', 'force' ]

    def validate(self):
        pass

    
class WavelengthSolutionPar(ParSet):
    def __init__(self, method=None, lamps=None, detection=None, numsearch=None, nfitpix=None,
                 IDpixels=None, IDwaves=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['method'] = 'arclines'
        options['method'] = WavelengthSolutionPar.valid_methods()
        dtypes['method'] = basestring
        descr['method'] = 'Method to use to fit the individual arc lines.  ' \
                          '\'fit\' is likely more accurate, but \'simple\' uses a polynomial ' \
                          'fit (to the log of a gaussian) and is fast and reliable.  ' \
                          '\'arclines\' uses the arclines python package.' \
                          'Options are: {0}'.format(', '.join(options['method']))

        # Force lamps to be a list
        if pars['lamps'] is not None and not isinstance(pars['lamps'], list):
            pars['lamps'] = [pars['lamps']]
        options['lamps'] = WavelengthSolutionPar.valid_lamps()
        dtypes['lamps'] = list
        descr['lamps'] = 'Name of one or more ions used for the wavelength calibration.  Use ' \
                         'None for no calibration.  ' \
                         'Options are: {0}'.format(', '.join(options['lamps']))

        defaults['detection'] = 6.0
        dtypes['detection'] = [int, float]
        descr['detection'] = 'Detection threshold for arc lines (in standard deviation)'

        defaults['numsearch'] = 20
        dtypes['numsearch'] = int
        descr['numsearch'] = 'Number of brightest arc lines to search for in preliminary ' \
                             'identification'

        defaults['nfitpix'] = 5
        dtypes['nfitpix'] = int
        descr['nfitpix'] = 'Number of pixels to fit when deriving the centroid of the arc ' \
                           'lines (an odd number is best)'

        dtypes['IDpixels'] = [int, list]
        descr['IDpixels'] = 'One or more pixels at which to manually identify a line'

        dtypes['IDwaves'] = [int, float, list]
        descr['IDwaves'] = 'Wavelengths of the manually identified lines'

        # Instantiate the parameter set
        super(WavelengthSolutionPar, self).__init__(list(pars.keys()),
                                                    values=list(pars.values()),
                                                    defaults=list(defaults.values()),
                                                    options=list(options.values()),
                                                    dtypes=list(dtypes.values()),
                                                    descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'method', 'lamps', 'detection', 'numsearch', 'nfitpix', 'IDpixels', 'IDwaves' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_methods():
        """
        Return the valid wavelength solution methods.
        """
        return [ 'simple', 'fit', 'arclines' ]

    @staticmethod
    def valid_lamps():
        """
        Return the valid lamp ions
        """
        return [ 'ArI', 'CdI', 'HgI', 'HeI', 'KrI', 'NeI', 'XeI', 'ZnI', 'ThAr' ]

    def validate(self):
        pass


class TraceSlitsPar(ParSet):
    """
    Parameters specific to PypIts slit tracing algorithm
    """
    def __init__(self, function=None, polyorder=None, medrep=None, number=None, trim=None,
                 maxgap=None, maxshift=None, pad=None, sigdetect=None, fracignore=None,
                 diffpolyorder=None, single=None, sobel_mode=None, pca=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        defaults['function'] = 'legendre'
        options['function'] = TraceSlitsPar.valid_functions()
        dtypes['function'] = basestring
        descr['function'] = 'Function use to trace the slit center.  ' \
                            'Options are: {0}'.format(', '.join(options['function']))

        defaults['polyorder'] = 3
        dtypes['polyorder'] = int
        descr['polyorder'] = 'Order of the function to use.'

        defaults['medrep'] = 0
        dtypes['medrep'] = int
        descr['medrep'] = 'Number of times to median smooth a trace image prior to analysis ' \
                          'for slit/order edges'

        # Force number to be an integer
        if values['number'] == 'auto':
            values['number'] = -1
        defaults['number'] = -1
        dtypes['number'] = int
        descr['number'] = 'Manually set the number of slits to identify (>=1). \'auto\' or -1 ' \
                          'will automatically identify the number of slits.'

        # Force trim to be a tuple
        if pars['trim'] is not None and not isinstance(pars['trim'], tuple):
            try:
                pars['trim'] = tuple(pars['trim'])
            except:
                raise TypeError('Could not convert provided trim to a tuple.')
        defaults['trim'] = (3,3)
        dtypes['trim'] = tuple
        descr['trim'] = 'How much to trim off each edge of each slit'

        dtypes['maxgap'] = int
        descr['maxgap'] = 'Maximum number of pixels to allow for the gap between slits.  Use ' \
                          'None if the neighbouring slits are far apart or of similar ' \
                          'illumination.'

        defaults['maxshift'] = 0.15
        dtypes['maxshift'] = [int, float]
        descr['maxshift'] = 'Maximum shift in trace crude'

        defaults['pad'] = 0
        dtypes['pad'] = int
        descr['pad'] = 'Integer number of pixels to consider beyond the slit edges.'

        defaults['sigdetect'] = 20.0
        dtypes['sigdetect'] = [int, float]
        descr['sigdetect'] = 'Sigma detection threshold for edge detection'
    
        defaults['fracignore'] = 0.01
        dtypes['fracignore'] = float
        descr['fracignore'] = 'If a slit spans less than this fraction over the spectral size ' \
                              'of the detector, it will be ignored (and reconstructed when/if ' \
                              'an \'order\' PCA analysis is performed).'

        defaults['diffpolyorder'] = 2
        dtypes['diffpolyorder'] = int
        descr['diffpolyorder'] = 'Order of the 2D function used to fit the 2d solution for the ' \
                                 'spatial size of all orders.'

        # TODO: Add a check for this?
        dtypes['single'] = list
        descr['single'] = 'Add a single, user-defined slit based on its location on each ' \
                          'detector.  Syntax is a list of values, 2 per detector, that define ' \
                          'the slit according to column values.  The second value (for the ' \
                          'right edge) must be greater than 0 to be applied.  LRISr example: ' \
                          'setting single = -1, -1, 7, 295 means the code will skip the ' \
                          'user-definition for the first detector but adds one for the second. ' \
                          ' None means no user-level slits defined.'

        defaults['sobel_mode'] = 'nearest'
        options['sobel_mode'] = TraceSlitsPar.valid_sobel_modes()
        dtypes['sobel_mode'] = basestring
        descr['sobel_mode'] = 'Mode for Sobel filtering.  Default is \'nearest\' but the ' \
                              'developers find \'constant\' works best for DEIMOS.'

        defaults['pca'] = PCAPar()
        dtypes['pca'] = [ ParSet, dict ]
        descr['pca'] = 'Parameters used for the PCA parameters in slit tracing'

        # Instantiate the parameter set
        super(TraceSlitsPar, self).__init__(list(pars.keys()),
                                            values=list(pars.values()),
                                            defaults=list(defaults.values()),
                                            options=list(options.values()),
                                            dtypes=list(dtypes.values()),
                                            descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'function', 'polyorder', 'medrep', 'number', 'trim', 'maxgap', 'maxshift',
                    'pad', 'sigdetect', 'fracignore', 'diffpolyorder', 'single', 'sobel_mode' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None

        pk = 'pca'
        kwargs[pk] = PCAPar.from_dict(cfg[pk]) if pk in k else None

        return cls(**kwargs)

    @staticmethod
    def valid_functions():
        """
        Return the list of valid functions to use for slit tracing.
        """
        return [ 'polynomial', 'legendre', 'chebyshev' ]

    @staticmethod
    def valid_sobel_modes():
        """Return the valid sobel modes."""
        return [ 'nearest', 'constant' ]

    def validate(self):
        if self.data['number'] == 0:
            raise ValueError('Number of slits must be -1 for automatic identification or '
                             'greater than 0')
        
class TraceTiltsPar(ParSet):
    """
    Parameters specific to PypIts slit tracing algorithm

    TODO: Changed to reflect wavetilts.py settings.  Was `yorder`
    previously `disporder`?  If so, I think I prefer the generality of
    `disporder`...
    """
    def __init__(self, idsonly=None, tracethresh=None, order=None, function=None, yorder=None,
                 func2D=None, method=None, params=None):


        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        defaults['idsonly'] = False
        dtypes['idsonly'] = bool
        descr['idsonly'] = 'Only use the arc lines that have an identified wavelength to trace ' \
                           'tilts'

        defaults['tracethresh'] = 1000.
        dtypes['tracethresh'] = [int, float]
        descr['tracethresh'] = 'TODO: X fill in the doc for this'

        defaults['order'] = 1
        dtypes['order'] = int
        descr['order'] = 'Order of the polynomial function to be used for the tilt of an ' \
                         'individual arc line.  Must be 1 for eschelle data (ARMED pipeline).'

        defaults['function'] = 'legendre'
        # TODO: Allowed values?
        dtypes['function'] = basestring
        descr['function'] = 'Type of function for arc line fits'

        defaults['yorder'] = 1
        dtypes['yorder'] = int
        descr['yorder'] = 'Order of the polynomial function to be used to fit the tilts ' \
                          'along the y direction.  TODO: Only used by ARMED pipeline?'

        defaults['func2D'] = 'legendre'
        # TODO: Allowed values?
        dtypes['func2D'] = basestring
        descr['func2D'] = 'Type of function for 2D fit'

        defaults['method'] = 'spline'
        options['method'] = TraceTiltsPar.valid_methods()
        dtypes['method'] = basestring
        descr['method'] = 'Method used to trace the tilt of the slit along an order.  ' \
                          'Options are: {0}'.format(', '.join(options['method']))

        # TODO: Need to add checks that check params against method
        defaults['params'] = [ 1, 1, 0 ]
        dtypes['params'] = [ int, list ]
        descr['params'] = 'Parameters to use for the provided method.  TODO: Need more explanation'

        # Instantiate the parameter set
        super(TraceTiltsPar, self).__init__(list(pars.keys()),
                                            values=list(pars.values()),
                                            defaults=list(defaults.values()),
                                            options=list(options.values()),
                                            dtypes=list(dtypes.values()),
                                            descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
#        parkeys = [ 'idsonly', 'method', 'params', 'order', 'disporder' ]
        parkeys = [ 'idsonly', 'tracethresh', 'order', 'function', 'yorder', 'func2D',
                    'method', 'params' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_methods():
        """
        Return the valid methods to use for tilt tracing.
        """
        return [ 'pca', 'spca', 'spline', 'interp', 'perp', 'zero' ]

    def validate(self):
        pass


class TraceObjectsPar(ParSet):
    """
    Parameters specific to PypIts object tracing algorithm
    """
    def __init__(self, function=None, order=None, find=None, nsmooth=None, xedge=None, method=None,
                 params=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['function'] = 'legendre'
        options['function'] = TraceObjectsPar.valid_functions()
        dtypes['function'] = basestring
        descr['function'] = 'Function to use to trace the object in each slit.  ' \
                            'Options are: {0}'.format(options['function'])

        defaults['order'] = 2
        dtypes['order'] = int
        descr['order'] = 'Order of the function to use to fit the object trace in each slit'

        defaults['find'] = 'standard'
        options['find'] = TraceObjectsPar.valid_detection_algorithms()
        dtypes['find'] = basestring
        descr['find'] = 'Algorithm to use for finding objects.' \
                        'Options are: {0}'.format(', '.join(options['find']))

        defaults['nsmooth'] = 3
        dtypes['nsmooth'] = [int, float]
        descr['nsmooth'] = 'Parameter for Gaussian smoothing when find=nminima.'

        defaults['xedge'] = 0.03
        dtypes['xedge'] = float
        descr['xedge'] = 'Ignore any objects within xedge of the edge of the slit'

        defaults['method'] = 'pca'
        options['method'] = TraceObjectsPar.valid_methods()
        dtypes['method'] = basestring
        descr['method'] = 'Method to use for tracing each object; only used with ARMED pipeline.' \
                          '  Options are: {0}'.format(', '.join(options['method']))

        defaults['params'] = [1, 0]
        dtypes['params'] = [int, list]
        descr['params'] = 'Parameters for the requested method.  For pca, params is a list ' \
                          'containing the order of the polynomials that should be used to fit ' \
                          'the object trace principal components. For example, params = 1, 0 ' \
                          'will fit 2 principal components, the first PC will be fit with a ' \
                          'first order polynomial, the second PC will be fit with a zeroth ' \
                          'order polynomial. TODO: What about the other methods?'

        # Instantiate the parameter set
        super(TraceObjectsPar, self).__init__(list(pars.keys()),
                                              values=list(pars.values()),
                                              defaults=list(defaults.values()),
                                              options=list(options.values()),
                                              dtypes=list(dtypes.values()),
                                              descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'function', 'order', 'find', 'nsmooth', 'xedge', 'method', 'params' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_functions():
        """
        Return the list of valid functions to use for object tracing.
        """
        return [ 'polynomial', 'legendre', 'chebyshev' ]

    @staticmethod
    def valid_detection_algorithms():
        """
        Return the list of valid algorithms for detecting objects.
        """
        return [ 'standard', 'nminima' ]

    @staticmethod
    def valid_methods():
        """
        Return the valid methods to use for tilt tracing.
        """
        return [ 'pca', 'spca', 'spline', 'interp', 'perp', 'zero' ]

    def validate(self):
        pass


class ExtractObjectsPar(ParSet):
    """
    Parameters specific to PypIts extraction of 1D object spectra
    """
    def __init__(self, pixelmap=None, pixelwidth=None, reuse=None, profile=None, maxnumber=None,
                 manual=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Check the manual input
        if manual is not None:
            if not isinstance(manual, (ParSet, dict, list)):
                raise TypeError('Manual extraction input must be a ParSet, dictionary, or list.')
            _manual = [manual] if isinstance(manual, (ParSet,dict)) else manual
            pars['manual'] = _manual

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        dtypes['pixelmap'] = basestring
        descr['pixelmap'] = 'If desired, a fits file can be specified (of the appropriate form)' \
                            'to specify the locations of the pixels on the detector (in ' \
                            'physical space).  TODO: Where is "appropriate form" specified?'

        defaults['pixelwidth'] = 2.5
        dtypes['pixelwidth'] = [int, float]
        descr['pixelwidth'] = 'The size of the extracted pixels (as an scaled number of Arc ' \
                              'FWHM), -1 will not resample'

        defaults['reuse'] = False
        dtypes['reuse'] = bool
        descr['reuse'] = 'If the extraction has previously been performed and saved, load the ' \
                         'previous result'

        defaults['profile'] = 'gaussian'
        options['profile'] = ExtractObjectsPar.valid_profiles()
        dtypes['profile'] = basestring
        descr['profile'] = 'Fitting function used to extract science data, only if the ' \
                           'extraction is 2D.  NOTE: options with suffix \'func\' fits a ' \
                           'function to the pixels whereas those without this suffix take into ' \
                           'account the integration of the function over the pixel (and is ' \
                           'closer to truth).   ' \
                           'Options are: {0}'.format(', '.join(options['profile']))

        dtypes['maxnumber'] = int
        descr['maxnumber'] = 'Maximum number of objects to extract in a science frame.  Use ' \
                             'None for no limit.'

        dtypes['manual'] = list
        descr['manual'] = 'List of manual extraction parameter sets'

        # Instantiate the parameter set
        super(ExtractObjectsPar, self).__init__(list(pars.keys()),
                                                values=list(pars.values()),
                                                defaults=list(defaults.values()),
                                                options=list(options.values()),
                                                dtypes=list(dtypes.values()),
                                                descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'pixelmap', 'pixelwidth', 'reuse', 'profile', 'maxnumber' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        kwargs['manual'] = _get_parset_list(cfg, 'manual', ManualExtractionPar)
        return cls(**kwargs)

    @staticmethod
    def valid_profiles():
        """
        Return the list of valid functions to use for object tracing.
        """
        return [ 'gaussian', 'gaussfunc', 'moffat', 'moffatfunc' ]

    def validate(self):
        pass

#-----------------------------------------------------------------------------
# Instrument ParSets

class DetectorPar(ParSet):
    def __init__(self, dataext=None, datasec=None, oscansec=None, dispaxis=None, xgap=None,
                 ygap=None, ysize=None, platescale=None, darkcurr=None, saturation=None,
                 nonlinear=None, numamplifiers=None, gain=None, ronoise=None, suffix=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['dataext'] = 0
        dtypes['dataext'] = int
        descr['dataext'] = 'Index of fits extension containing data'

        # TODO: Allow for None, such that the entire image is the data
        # section
        defaults['datasec'] = 'DATASEC'
        dtypes['datasec'] = basestring
        descr['datasec'] = 'Either the data sections or the header keyword where the valid ' \
                           'data sections can be obtained. If defined explicitly should have ' \
                           'the format of a numpy array slice'

        # TODO: Allow for None, such that there is no overscan region
        defaults['oscansec'] = 'BIASSEC'
        dtypes['oscansec'] = basestring
        descr['oscansec'] = 'Either the overscan section or the header keyword where the valid ' \
                            'data sections can be obtained. If defined explicitly should have ' \
                            'the format of a numpy array slice'

        # TODO: Should this be detector-specific, or camera-specific?
        defaults['dispaxis'] = 0
        options['dispaxis'] = [ 0, 1 ]
        dtypes['dispaxis'] = int
        descr['dispaxis'] = 'Spectra are dispersed along this axis (0 for row, 1 for column)'

        defaults['xgap'] = 0.0
        dtypes['xgap'] = [int, float]
        descr['xgap'] = 'Gap between the square detector pixels (expressed as a fraction of the ' \
                        'x pixel size -- x is predominantly the dispersion axis)'

        defaults['ygap'] = 0.0
        dtypes['ygap'] = [int, float]
        descr['ygap'] = 'Gap between the square detector pixels (expressed as a fraction of the ' \
                        'y pixel size -- x is predominantly the dispersion axis)'

        defaults['ysize'] = 1.0
        dtypes['ysize'] = [int, float]
        descr['ysize'] = 'The size of a pixel in the y-direction as a multiple of the x pixel ' \
                         'size (i.e. xsize = 1.0 -- x is predominantly the dispersion axis)'

        defaults['platescale'] = 0.135
        dtypes['platescale'] = [int, float]
        descr['platescale'] = 'arcsec per pixel in the spatial dimension for an unbinned pixel'

        defaults['darkcurr'] = 0.0
        dtypes['darkcurr'] = [int, float]
        descr['darkcurr'] = 'Dark current (e-/hour)'

        defaults['saturation'] = 65535.0
        dtypes['saturation'] = [ int, float ]
        descr['saturation'] = 'The detector saturation level'

        defaults['nonlinear'] = 0.86
        dtypes['nonlinear'] = [ int, float ]
        descr['nonlinear'] = 'Percentage of detector range which is linear (i.e. everything ' \
                             'above nonlinear*saturation will be flagged as saturated)'

        defaults['numamplifiers'] = 1
        dtypes['numamplifiers'] = int
        descr['numamplifiers'] = 'Number of amplifiers'

        defaults['gain'] = 1.0
        dtypes['gain'] = [ int, float, list ]
        descr['gain'] = 'Inverse gain (e-/ADU). A list should be provided if a detector ' \
                        'contains more than one amplifier.'

        defaults['ronoise'] = 4.0
        dtypes['ronoise'] = [ int, float, list ]
        descr['ronoise'] = 'Read-out noise (e-). A list should be provided if a detector ' \
                           'contains more than one amplifier.'

        # TODO: Allow this to be None?
        defaults['suffix'] = ''
        dtypes['suffix'] = basestring
        descr['suffix'] = 'Suffix to be appended to all saved calibration and extraction frames.'

        # Instantiate the parameter set
        super(DetectorPar, self).__init__(list(pars.keys()),
                                          values=list(pars.values()),
                                          defaults=list(defaults.values()),
                                          options=list(options.values()),
                                          dtypes=list(dtypes.values()),
                                          descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'dataext', 'datasec', 'oscansec', 'dispaxis', 'xgap', 'ygap', 'ysize',
                    'platescale', 'darkcurr', 'saturation', 'nonlinear', 'numamplifiers', 'gain',
                    'ronoise', 'suffix' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['numamplifiers'] > 1:
            if not isinstance(self.data['gain'], list) \
                    or len(self.data['gain']) != self.data['numamplifiers']:
                raise ValueError('Provided gain is not a list of the correct length.')
            if not isinstance(self.data['ronoise'], list) \
                    or len(self.data['ronoise']) != self.data['numamplifiers']:
                raise ValueError('Provided read-out noise is not a list of the correct length.')

            for k in range(self.data['numamplifiers']):
                if not isinstance(self.data['gain'][k], (int, float)):
                    TypeError('Gain values must be a integer or floating-point number.')
                if not isinstance(self.data['ronoise'][k], (int, float)):
                    TypeError('Read-out noise values must be a integer or floating-point number.')


class InstrumentPar(ParSet):
    def __init__(self, telescope=None, longitude=None, latitude=None, elevation=None, camera=None,
                 minexp=None, detector=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Check the detector input
        if detector is not None:
            if not isinstance(detector, (ParSet, dict, list)):
                raise TypeError('Detector input must be a ParSet, dictionary, or list.')
            _detector = [detector] if isinstance(detector, (ParSet,dict)) else detector
            pars['detector'] = _detector

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['telescope'] = 'KECK'
        options['telescope'] = InstrumentPar.valid_telescopes()
        dtypes['telescope'] = basestring
        descr['telescope'] = 'Name of the telescope used to obtain the observations.  ' \
                        'Options are: {0}'.format(', '.join(options['telescope']))
        
        defaults['longitude'] = 155.47833
        dtypes['longitude'] = [int, float]
        descr['longitude'] = 'Longitude of the telescope on Earth in degrees.'

        defaults['latitude'] = 19.82833
        dtypes['latitude'] = [int, float]
        descr['latitude'] = 'Latitude of the telescope on Earth in degrees.'

        defaults['elevation'] = 4160.0
        dtypes['elevation'] = [int, float]
        descr['elevation'] = 'Elevation of the telescope in m'

        defaults['camera'] = 'LRISb'
        options['camera'] = InstrumentPar.valid_cameras()
        dtypes['camera'] = basestring
        descr['camera'] = 'Name of the camera.  ' \
                        'Options are: {0}'.format(', '.join(options['camera']))
        
        defaults['minexp'] = 1.0
        dtypes['minexp'] = [int, float]
        descr['minexp'] = 'Minimum exposure time in seconds'

        defaults['detector'] = [ DetectorPar() ]
        dtypes['detector'] = list
        descr['detector'] = 'List of detectors'

        # Instantiate the parameter set
        super(InstrumentPar, self).__init__(list(pars.keys()),
                                            values=list(pars.values()),
                                            defaults=list(defaults.values()),
                                            options=list(options.values()),
                                            dtypes=list(dtypes.values()),
                                            descr=list(descr.values()))

        # Check the parameters match the method requirements
        self.validate()


    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'telescope', 'longitude', 'latitude', 'elevation', 'camera', 'minexp' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        kwargs['detector'] = _get_parset_list(cfg, 'detector', DetectorPar)
        return cls(**kwargs)

    @staticmethod
    def valid_telescopes():
        """
        Return the valid telescopes.
        """
        return [ 'KECK', 'SHANE', 'WHT', 'APF', 'TNG' ]

    @staticmethod
    def valid_cameras():
        """
        Return the list of valid cameras.

        TODO: This should probably be telescope dependent...
        """
        return [ 'LEVY', 'DEIMOS', 'HIRES', 'LRISb', 'LRISr', 'NIRSPEC', 'KASTb', 'KASTr',
                 'DOLORES', 'ISISb' ]

    def validate(self):
        if self.data['detector'] is None or len(self.data['detector']) < 1:
            raise ValueError('Must defined at least one detector!')


class FrameFitsPar(ParSet):
    def __init__(self, timeunit=None, headext=None, lamps=None, keydef=None, keycheck=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['timeunit'] = 'mjd'
        options['timeunit'] = FrameFitsPar.valid_time_units()
        dtypes['timeunit'] = basestring
        descr['timeunit'] = 'The unit of time keyword.  ' \
                            'Options are: {0}'.format(', '.join(options['timeunit']))

        defaults['headext'] = 0
        dtypes['headext'] = [int, list]
        descr['headext'] = '(List of) Extension(s) with data to read (0-indexed)'

        dtypes['lamps'] = [basestring, list]
        descr['lamps'] = 'One or more lamp names'

        # Force all the keyword definitions to provide the associated
        # extension
        if pars['keydef'] is not None:
            for k in pars['keydef'].keys():
                if not isinstance(pars['keydef'][k], list):
                    pars['keydef'][k] = [0, pars['keydef'][k]]
        dtypes['keydef'] = dict
        descr['keydef'] = 'Dictionary with the definitions of keywords to used from the ' \
                          'header.  Variable names must be unique and define a unique header ' \
                          'keyword.  If a single string value, the keyword is expected to be ' \
                          'in the primary (extension 0) header.  If defined as a 2-element ' \
                          'list, the first element is the extension with the approipriate ' \
                          'header keyword (0-indexed) and the second is the keyword name.'

        dtypes['keycheck'] = dict
        descr['keycheck'] = 'Dictionary with checks to perform on the keyword values for all ' \
                            'fits files.  The dictionary keyword must be one of the defined ' \
                            'header keywords in \'keydef\'.  Single values imply an equality ' \
                            'check.  A list with two values implies a lower and upper range.  ' \
                            'To set *only* a lower or upper limit, set the unconstrained limit ' \
                            'to None.  All limits are exclusive.  I.e, to get a value >0, set ' \
                            '\'0, None\'; to get <30, set \'None, 30\'.'

        # Instantiate the parameter set
        super(FrameFitsPar, self).__init__(list(pars.keys()),
                                           values=list(pars.values()),
                                           defaults=list(defaults.values()),
                                           options=list(options.values()),
                                           dtypes=list(dtypes.values()),
                                           descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'timeunit', 'headext', 'lamps', 'keydef', 'keycheck' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_time_units():
        """
        Return the valid time units.
        """
        return [ 'h', 'm', 's' ] + list(Time.FORMATS.keys())

    def validate(self):
        """Validate the parameter set."""
        # Check the keyword definitions make sense
        if self.data['keydef'] is not None:
            for k in self.data['keydef'].keys():
                if not isinstance(self.data['keydef'][k], list):
                    raise TypeError('The definition of {0} is not a list!'.format(k))
                if len(self.data['keydef'][k]) != 2:
                    raise ValueError('The list definition of {0} is too long!'.format(k))
                if not isinstance(self.data['keydef'][k][0], int):
                    raise TypeError('The extension for keyword {0} must be an int!'.format(k))
                if not isinstance(self.data['keydef'][k][1], basestring):
                    raise TypeError('The header keyword for {0} must be a string!'.format(k))

        # Only keyword definitions were provided
        if self.data['keycheck'] is None:
            return

        # Cannot provide keyword checks without the keyword definitions
        if self.data['keydef'] is None and self.data['keycheck'] is not None:
            raise ValueError('To apply keyword checks, must first define a set of keywords!')

        # Check all keyword checks have an associated keyword definition
        for k in self.data['keycheck'].keys():
            if k not in self.data['keydef'].keys():
                raise KeyError('{0} does not have an associated keyword definition.'.format(k))


class FrameIDPar(ParSet):
    def __init__(self, fitspar=None, frametype=None, canbe=None, keycheck=None, match=None):

        # Save the pointer to the FrameFitsPar
        self.fitspar = FrameFitsPar() if fitspar is None else fitspar
        if not isinstance(self.fitspar, (ParSet,dict)):
            raise TypeError('fitspar must be defined in FrameIDPar must be a ParSet or dict.')

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[2:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['frametype'] = 'bias'
        options['frametype'] = FrameGroupPar.valid_frame_types()
        dtypes['frametype'] = basestring
        descr['frametype'] = 'Frame type.  ' \
                             'Options are: {0}'.format(', '.join(options['frametype']))

        dtypes['canbe'] = basestring
        descr['canbe'] = 'Frame types, **other than this type**, that can also be used as ' \
                         'this type.  E.g., this frametype is \'pixelflat\' but can also be '\
                         'used as a \'trace\' frame.'

        dtypes['keycheck'] = dict
        descr['keycheck'] = 'Dictionary with checks to perform on the keyword values to select ' \
                            'frames of this type.  The dictionary keyword must be one of the ' \
                            'defined header keywords in the associated \`fitspar\` parameter ' \
                            'set.  Single values imply an equality check.  A list with two ' \
                            'values implies a lower and upper range.  To set *only* a lower or ' \
                            'upper limit, set the unconstrained limit to None.  All limits are ' \
                            'exclusive.  I.e, to get a value >0, set \'0, None\'; to get <30, ' \
                            'set \'None, 30\'.'

        dtypes['match'] = dict
        descr['match'] = 'Dictionary with properties used to isolate frames of this type.  ' \
                         'TODO: Give examples of how to use this'
        
        # Instantiate the parameter set
        super(FrameIDPar, self).__init__(list(pars.keys()),
                                         values=list(pars.values()),
                                         defaults=list(defaults.values()),
                                         options=list(options.values()),
                                         dtypes=list(dtypes.values()),
                                         descr=list(descr.values()))

        self.validate()

    @classmethod
    def from_dict(cls, fitspar, frametype, cfg):
        k = cfg.keys()
        parkeys = [ 'canbe', 'keycheck', 'match' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(fitspar=fitspar, frametype=frametype, **kwargs)

    def validate(self):
        """Validate the parameter set."""
        # No checks are provided
        if self.data['keycheck'] is None:
            return

        # Cannot provide keyword checks without the keyword definitions
        if self.fitspar['keydef'] is None and self.data['keycheck'] is not None:
            raise ValueError('To apply keyword checks, must define fits keywords!')

        # Check all keyword checks have an associated keyword definition
        for k in self.data['keycheck'].keys():
            if k == 'lampstat':
                # lampstat is a special keyword that doesn't
                # (necessarily) have a directly associated header
                # keyword
                continue
            if k not in self.fitspar['keydef'].keys():
                raise KeyError('{0} does not have an associated keyword definition.'.format(k))


#-----------------------------------------------------------------------------
# Parameters superset

class PypitPar(ParSet):
    """
    The superset of all parameters used by PypIt.

    Users will likely always want to use the :func:`from_cfg_file`
    method to instantiate the parameter set, instead of using this
    instantiation function.

    .. todo::

        - Is there a better way we can identify and group frames?  Do we
          need to carry around the *id and *group parameters during the
          entire pypit run?
        - Should the FrameIDPar groups become part of InstrumentPar?
    """
    def __init__(self, run=None, rdx=None, instrument=None, fits=None, biasid=None,
                 pixelflatid=None, arcid=None, pinholeid=None, traceid=None, standardid=None,
                 scienceid=None, biasgroup=None, pixelflatgroup=None, arcgroup=None,
                 pinholegroup=None, tracegroup=None, standardgroup=None, sciencegroup=None,
                 wavelengths=None, slits=None, tilts=None, objects=None, extract=None):

        # Components set internally by the code and not by the user

        # TODO: Not sure we want these here
        try:
            self.calling_program = __file__     # Name of the calling program
        except NameError:
            self.calling_program = None

        # PypIt root directory
        # TODO: This should go in a different module
        self.pypit_root = _pypit_root_directory()
        
        # TODO: Not sure this is needed
        self.user_par = None                    # File with the user-defined parameters

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['run'] = RunPar()
        dtypes['run'] = [ ParSet, dict ]
        descr['run'] = 'PypIt execution options.'

        defaults['rdx'] = ReducePar()
        dtypes['rdx'] = [ ParSet, dict ]
        descr['rdx'] = 'PypIt reduction rules.'

        # TODO: Should there be meaningful default instrument, fits, and id
        # parameter sets?  The full parameter set is meaningless without
        # them, and the current defaults just set the keywords for each
        # sub-parameter set.
        defaults['instrument'] = InstrumentPar()
        dtypes['instrument'] = [ ParSet, dict ]
        descr['instrument'] = 'PypIt instrument parameters.'

        defaults['fits'] = FrameFitsPar()
        dtypes['fits'] = [ ParSet, dict ]
        descr['fits'] = 'The fits file parameters and checks general to all data from the ' \
                        'instrument to be reduced'

        defaults['biasid'] = FrameIDPar(fitspar=defaults['fits'], frametype='bias')
        dtypes['biasid'] = [ ParSet, dict ]
        descr['biasid'] = 'The identification rules and checks for bias frames'

        defaults['pixelflatid'] = FrameIDPar(fitspar=defaults['fits'], frametype='pixelflat')
        dtypes['pixelflatid'] = [ ParSet, dict ]
        descr['pixelflatid'] = 'The identification rules and checks for pixel-flat frames'

        defaults['arcid'] = FrameIDPar(fitspar=defaults['fits'], frametype='arc')
        dtypes['arcid'] = [ ParSet, dict ]
        descr['arcid'] = 'The identification rules and checks for arc frames'

        defaults['pinholeid'] = FrameIDPar(fitspar=defaults['fits'], frametype='pinhole')
        dtypes['pinholeid'] = [ ParSet, dict ]
        descr['pinholeid'] = 'The identification rules and checks for pin-hole frames'

        defaults['traceid'] = FrameIDPar(fitspar=defaults['fits'], frametype='trace')
        dtypes['traceid'] = [ ParSet, dict ]
        descr['traceid'] = 'The identification rules and checks for trace frames'

        defaults['standardid'] = FrameIDPar(fitspar=defaults['fits'], frametype='standard')
        dtypes['standardid'] = [ ParSet, dict ]
        descr['standardid'] = 'The identification rules and checks for standard frames'

        defaults['scienceid'] = FrameIDPar(fitspar=defaults['fits'], frametype='science')
        dtypes['scienceid'] = [ ParSet, dict ]
        descr['scienceid'] = 'The identification rules and checks for science frames'

        defaults['biasgroup'] = FrameGroupPar(frametype='bias')
        dtypes['biasgroup'] = [ ParSet, dict ]
        descr['biasgroup'] = 'The frames and combination rules for the bias correction'

        defaults['pixelflatgroup'] = FrameGroupPar(frametype='pixelflat')
        dtypes['pixelflatgroup'] = [ ParSet, dict ]
        descr['pixelflatgroup'] = 'The frames and combination rules for the field flattening'

        defaults['arcgroup'] = FrameGroupPar(frametype='arc')
        dtypes['arcgroup'] = [ ParSet, dict ]
        descr['arcgroup'] = 'The frames and combination rules for the wavelength calibration'

        defaults['pinholegroup'] = FrameGroupPar(frametype='pinhole')
        dtypes['pinholegroup'] = [ ParSet, dict ]
        descr['pinholegroup'] = 'The frames and combination rules for tracing the slit centroid'

        defaults['tracegroup'] = FrameGroupPar(frametype='trace')
        dtypes['tracegroup'] = [ ParSet, dict ]
        descr['tracegroup'] = 'The frames and combination rules for tracing the slit edges'

        defaults['standardgroup'] = FrameGroupPar(frametype='standard')
        dtypes['standardgroup'] = [ ParSet, dict ]
        descr['standardgroup'] = 'The frames and combination rules for the spectrophotometric ' \
                                 'standard observations'

        defaults['sciencegroup'] = FrameGroupPar(frametype='science')
        dtypes['sciencegroup'] = [ ParSet, dict ]
        descr['sciencegroup'] = 'The frames and combination rules for the science observations'

        defaults['wavelengths'] = WavelengthSolutionPar()
        dtypes['wavelengths'] = [ ParSet, dict ]
        descr['wavelengths'] = 'Parameters used to derive the wavelength solution'

        defaults['slits'] = TraceSlitsPar()
        dtypes['slits'] = [ ParSet, dict ]
        descr['slits'] = 'Define how the slits should be traced using the trace ?PINHOLE? frames'

        defaults['tilts'] = TraceTiltsPar()
        dtypes['tilts'] = [ ParSet, dict ]
        descr['tilts'] = 'Define how to tract the slit tilts using the trace frames'

        defaults['objects'] = TraceObjectsPar()
        dtypes['objects'] = [ ParSet, dict ]
        descr['objects'] = 'Define how to tract the slit tilts using the trace frames'

        defaults['extract'] = ExtractObjectsPar()
        dtypes['extract'] = [ ParSet, dict ]
        descr['extract'] = 'Define how to extract 1D object spectra'

        # Instantiate the parameter set
        super(PypitPar, self).__init__(list(pars.keys()),
                                       values=list(pars.values()),
                                       defaults=list(defaults.values()),
                                       dtypes=list(dtypes.values()),
                                       descr=list(descr.values()))

        self.validate()

    @classmethod
    def from_cfg_file(cls, cfg_file=None, merge_with=None, expand_spectrograph=True,
                      evaluate=True):
        """
        Construct the parameter set using a configuration file.

        Note that::

            default = PypitPar()
            nofile = PypitPar.from_cfg_file()
            assert default.data == nofile.data, 'This should always pass.'

        Args:
            cfg_file (:obj:`str`, optional):
                The name of the configuration file that defines the
                default parameters.  This can be used if have a pypit
                config file from a previous run that was constructed and
                output by pypit.  This has to contain the full set of
                parameters, not just the subset you want to change.  For
                the latter, use :arg:`merge_with` to provide one or more
                config files to merge with the defaults to construct the
                full parameter set.
            merge_with (:obj:`str`, :obj:`list`, optional):
                One or more config files with the modifications to
                either default parameters (:arg:`cfg_file` is None) or
                the parameters provided by :arg:`cfg_file`.  The
                modifications are performed in series so the list order
                of the config files is important.
            expand_spectrograph (:obj:`bool`, optional):
                Use the `cfg['rdx']['spectrograph']` keyword to select
                and expand the instrument parameters using a
                configuration file called::

                   '{0}_spectrograph.cfg'.format(cfg['rdx']['spectrograph'])

                The name of the spectrograph is defined by the merged
                sequence of config files, following that merging
                precendence (see above).  A ValueError is raised if the
                merged config values have::
                
                    cfg['rdx']['spectrograph'] == 'None' 

                Once the spectrograph configuration is read, the merging
                sequence is as follows: (1) the default configuration
                (set by the :arg:`cfg_file` argument), (2) the default
                spectrograph parameters read by this keyword selection,
                and then (3) the modifications set by the merging
                sequence.  This allows the user to select the default
                spectrograph configuration and then alter any of the
                parameters defined in either the reduction or
                spectrograph sets in one config file.
            evaluate (:obj:`bool`, optional):
                Evaluate the values in the config object before
                assigning them in the subsequent parameter sets.  The
                parameters in the config file are *always* read as
                strings, so this should almost always be true; however,
                see the warning below.
                
        .. warning::

            When :arg:`evaluate` is true, the function runs `eval()` on
            all the entries in the `ConfigObj` dictionary, done using
            :func:`_recursive_dict_evaluate`.  This has the potential to
            go haywire if the name of a parameter unintentionally
            happens to be identical to an imported or system-level
            function.  Of course, this can be useful by allowing one to
            define the function to use as a parameter, but it also means
            one has to be careful with the values that the parameters
            should be allowed to have.  The current way around this is
            to provide a list of strings that should be ignored during
            the evaluation, done using :func:`_eval_ignore`.

        Returns:
            :class:`pypit.par.core.PypitPar`: The instance of the
            parameter set.

        Raises:
            ValueError: Raised if the spectrograph keyword is 'None' and
                `expand_spectrograph=True`.

        """
        # Get the base parameters in a ConfigObj instance
        cfg = ConfigObj(PypitPar().to_config(None, just_lines=True)
                            if cfg_file is None else cfg_file)

        # Get the list of other configuration parameters to merge it with
        _merge_with = [] if merge_with is None else \
                        ([merge_with] if isinstance(merge_with, basestring) else merge_with)
        merge_cfg = ConfigObj()
        for f in _merge_with:
            merge_cfg.merge(ConfigObj(f))

        # Use the keyword set for the spectrograph to grab the
        # spectrograph configuration
        spec_cfg = ConfigObj()
        if expand_spectrograph:
            try:
                spectrograph = merge_cfg['rdx']['spectrograph']
            except:
                spectrograph = 'None'
            if spectrograph == 'None':
                spectrograph = cfg['rdx']['spectrograph']
            if spectrograph == 'None':
                raise ValueError('Spectrograph is undefined!')
            spectrograph_cfg_file = ReducePar.spectrograph_config_file(spectrograph, verbose=True)
            spec_cfg = ConfigObj(spectrograph_cfg_file)

        # The merge order is default, spectrograph, merge.  The merge
        # will be successful if either spec_cfg or merge_cfg are empty
        # ConfigObj instances.
        cfg.merge(spec_cfg)
        cfg.merge(merge_cfg)

#        cfg = _recursive_dict_unicode2str(cfg)

        # Evaluate the strings if requested
        if evaluate:
            cfg = _recursive_dict_evaluate(cfg)
        
        # Instantiate the object based on the configuration dictionary
        return cls.from_dict(cfg)

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        kwargs = {}

        pk = 'run'
        kwargs[pk] = RunPar.from_dict(cfg[pk]) if pk in k else None

        pk = 'rdx'
        kwargs[pk] = ReducePar.from_dict(cfg[pk]) if pk in k else None

        pk = 'instrument'
        kwargs[pk] = InstrumentPar.from_dict(cfg[pk]) if pk in k else None

        pk = 'fits'
        kwargs[pk] = FrameFitsPar.from_dict(cfg[pk]) if pk in k else None

        pk = 'biasid'
        kwargs[pk] = FrameIDPar.from_dict(kwargs['fits'], 'bias', cfg[pk]) if pk in k else None

        pk = 'pixelflatid'
        kwargs[pk] = FrameIDPar.from_dict(kwargs['fits'], 'pixelflat', cfg[pk]) if pk in k else None

        pk = 'arcid'
        kwargs[pk] = FrameIDPar.from_dict(kwargs['fits'], 'arc', cfg[pk]) if pk in k else None

        pk = 'pinholeid'
        kwargs[pk] = FrameIDPar.from_dict(kwargs['fits'], 'pinhole', cfg[pk]) if pk in k else None

        pk = 'traceid'
        kwargs[pk] = FrameIDPar.from_dict(kwargs['fits'], 'trace', cfg[pk]) if pk in k else None

        pk = 'standardid'
        kwargs[pk] = FrameIDPar.from_dict(kwargs['fits'], 'standard', cfg[pk]) if pk in k else None

        pk = 'scienceid'
        kwargs[pk] = FrameIDPar.from_dict(kwargs['fits'], 'science', cfg[pk]) if pk in k else None

        pk = 'biasgroup'
        kwargs[pk] = FrameGroupPar.from_dict('bias', cfg[pk]) if pk in k else None

        pk = 'pixelflatgroup'
        kwargs[pk] = FrameGroupPar.from_dict('pixelflat', cfg[pk]) if pk in k else None

        pk = 'arcgroup'
        kwargs[pk] = FrameGroupPar.from_dict('arc', cfg[pk]) if pk in k else None

        pk = 'pinholegroup'
        kwargs[pk] = FrameGroupPar.from_dict('pinhole', cfg[pk]) if pk in k else None

        pk = 'tracegroup'
        kwargs[pk] = FrameGroupPar.from_dict('trace', cfg[pk]) if pk in k else None

        pk = 'standardgroup'
        kwargs[pk] = FrameGroupPar.from_dict('standard', cfg[pk]) if pk in k else None

        pk = 'sciencegroup'
        kwargs[pk] = FrameGroupPar.from_dict('science', cfg[pk]) if pk in k else None

        pk = 'wavelengths'
        kwargs[pk] = WavelengthSolutionPar.from_dict(cfg[pk]) if pk in k else None

        pk = 'slits'
        kwargs[pk] = TraceSlitsPar.from_dict(cfg[pk]) if pk in k else None
        
        pk = 'tilts'
        kwargs[pk] = TraceTiltsPar.from_dict(cfg[pk]) if pk in k else None

        pk = 'objects'
        kwargs[pk] = TraceObjectsPar.from_dict(cfg[pk]) if pk in k else None

        pk = 'extract'
        kwargs[pk] = ExtractObjectsPar.from_dict(cfg[pk]) if pk in k else None

        return cls(**kwargs)

    # TODO: Perform extensive checking that the parameters are valid for
    # a full run of PYPIT.  May not be necessary because validate will
    # be called for all the sub parameter sets, but this can do higher
    # level checks, if necessary.
    def validate(self):
        pass

