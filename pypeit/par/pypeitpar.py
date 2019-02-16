# encoding: utf-8
"""
Defines parameter sets used to set the behavior for core pypeit
functionality.

For more details on the full parameter hierarchy and a tabulated
description of the keywords in each parameter set, see :ref:`pypeitpar`.

For examples of how to change the parameters for a run of pypeit using
the pypeit input file, see :ref:`pypeit_file`.

**New Parameters**:

To add a new parameter, let's call it `foo`, to any of the provided
parameter sets:

    - Add ``foo=None`` to the ``__init__`` method of the relevant
      parameter set.  E.g.::
        
        def __init__(self, existing_par=None, foo=None):

    - Add any default value (the default value is ``None`` unless you set
      it), options list, data type, and description to the body of the
      ``__init__`` method.  E.g.::

        defaults['foo'] = 'bar'
        options['foo'] = [ 'bar', 'boo', 'fighters' ]
        dtypes['foo'] = str
        descr['foo'] = 'foo? who you callin a foo!  ' \
                       'Options are: {0}'.format(', '.join(options['foo']))

    - Add the parameter to the ``from_dict`` method:
    
        - If the parameter is something that does not require
          instantiation, add the keyword to the ``parkeys`` list in the
          ``from_dict`` method.  E.g.::

            parkeys = [ 'existing_par', 'foo' ]
            kwargs = {}
            for pk in parkeys:
                kwargs[pk] = cfg[pk] if pk in k else None

        - If the parameter is another ParSet or requires instantiation,
          provide the instantiation.  For example, see how the
          :class:`ProcessImagesPar` parameter set is defined in the
          :class:`FrameGroupPar` class.  E.g.::

            pk = 'foo'
            kwargs[pk] = FooPar.from_dict(cfg[pk]) if pk in k else None

**New Parameter Sets:**

To add an entirely new parameter set, use one of the existing parameter
sets as a template, then add the parameter set to :class:`PypeItPar`,
assuming you want it to be accessed throughout the code.

----
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import glob
import warnings
from pkg_resources import resource_filename
import inspect

from collections import OrderedDict

import numpy

from configobj import ConfigObj
from astropy.time import Time

from pypeit.par.parset import ParSet
from pypeit.par import util

# Needs this to determine the valid spectrographs TODO: This causes a
# circular import.  Spectrograph specific parameter sets and where they
# go needs to be rethought.
#from ..spectrographs.util import valid_spectrographs

#-----------------------------------------------------------------------------
# Reduction ParSets

# TODO: Create child classes for each allowed frame type?  E.g.:
#
# class BiasPar(FrameGroupPar):
#    def __init__(self, useframe=None, number=None, overscan=None, combine=None, lacosmic=None):
#        # Set frame-specific defaults
#        _number = 5 if number is None else number
#        super(BiasPar, self).__init(frametype='bias', useframe=useframe, number=_number,
#                                    overscan=overscan, combine=combine, lacosmic=lacosmic)

class FrameGroupPar(ParSet):
    """
    An abstracted group of parameters that defines how specific types of
    frames should be grouped and combined.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """
    def __init__(self, frametype=None, useframe=None, number=None, exprng=None, process=None):
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
#        defaults['frametype'] = 'bias'
        defaults['frametype'] = frametype       # This is a kludge
        options['frametype'] = FrameGroupPar.valid_frame_types()
        dtypes['frametype'] = str
        descr['frametype'] = 'Frame type.  ' \
                             'Options are: {0}'.format(', '.join(options['frametype']))

        dtypes['useframe'] = str
        descr['useframe'] = 'A master calibrations file to use if it exists.'

        defaults['number'] = 0
        dtypes['number'] = int
        descr['number'] = 'Used in matching calibration frames to science frames.  This sets ' \
                          'the number of frames to use of this type'

        defaults['exprng'] = [None, None]
        dtypes['exprng'] = list
        descr['exprng'] = 'Used in identifying frames of this type.  This sets the minimum ' \
                          'and maximum allowed exposure times.  There must be two items in ' \
                          'the list.  Use None to indicate no limit; i.e., to select exposures ' \
                          'with any time greater than 30 sec, use exprng = [30, None].'

        defaults['process'] = ProcessImagesPar()
        dtypes['process'] = [ ParSet, dict ]
        descr['process'] = 'Parameters used for basic image processing'

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
        parkeys = [ 'useframe', 'number', 'exprng' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        pk = 'process'
        kwargs[pk] = ProcessImagesPar.from_dict(cfg[pk]) if pk in k else None
        return cls(frametype=frametype, **kwargs)

    @staticmethod
    def valid_frame_types():
        """
        Return the list of valid frame types.
        """
        return [ 'bias', 'dark', 'pixelflat', 'arc', 'pinhole', 'trace', 'standard', 'science',
                 'all' ]

    def validate(self):
        if self.data['useframe'] is None:
            self.default['useframe'] = self.data['frametype']
            self.data['useframe'] = self.data['frametype']
        if len(self.data['exprng']) != 2:
            raise ValueError('exprng must be a list with two items.')


class ProcessImagesPar(ParSet):
    """
    The parameters needed to perform basic image processing.

    These parameters are primarily used by
    :class:`pypeit.processimages.ProcessImages`, the base class of many
    of the pypeit objects.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """
    def __init__(self, overscan=None, overscan_par=None, match=None, combine=None, satpix=None,
                 sigrej=None, n_lohi=None, sig_lohi=None, replace=None, lamaxiter=None, grow=None,
                 rmcompact=None, sigclip=None, sigfrac=None, objlim=None):

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
        defaults['overscan'] = 'savgol'
        options['overscan'] = ProcessImagesPar.valid_overscan()
        dtypes['overscan'] = str
        descr['overscan'] = 'Method used to fit the overscan.  ' \
                            'Options are: {0}'.format(', '.join(options['overscan']))
        
        defaults['overscan_par'] = [5, 65]
        dtypes['overscan_par'] = [int, list]
        descr['overscan_par'] = 'Parameters for the overscan subtraction.  For ' \
                                '\'polynomial\', set overcan_par = order, number of pixels, ' \
                                'number of repeats ; for \'savgol\', set overscan_par = ' \
                                'order, window size ; for \'median\', set overscan_par = ' \
                                'None or omit the keyword.'

        defaults['match'] = -1
        dtypes['match'] = [int, float]
        descr['match'] = '(Deprecate?) Match frames with pixel counts that are within N-sigma ' \
                         'of one another, where match=N below.  If N < 0, nothing is matched.'

        defaults['combine'] = 'weightmean'
        options['combine'] = ProcessImagesPar.valid_combine_methods()
        dtypes['combine'] = str
        descr['combine'] = 'Method used to combine frames.  Options are: {0}'.format(
                                       ', '.join(options['combine']))

        defaults['satpix'] = 'reject'
        options['satpix'] = ProcessImagesPar.valid_saturation_handling()
        dtypes['satpix'] = str
        descr['satpix'] = 'Handling of saturated pixels.  Options are: {0}'.format(
                                       ', '.join(options['satpix']))

        defaults['sigrej'] = 20.0
        dtypes['sigrej'] = [int, float]
        descr['sigrej'] = 'Sigma level to reject cosmic rays (<= 0.0 means no CR removal)'

        defaults['n_lohi'] = [0, 0]
        dtypes['n_lohi'] = list
        descr['n_lohi'] = 'Number of pixels to reject at the lowest and highest ends of the ' \
                          'distribution; i.e., n_lohi = low, high.  Use None for no limit.'

        defaults['sig_lohi'] = [3.0, 3.0]
        dtypes['sig_lohi'] = list
        descr['sig_lohi'] = 'Sigma-clipping level at the low and high ends of the distribution; ' \
                            'i.e., sig_lohi = low, high.  Use None for no limit.'

        defaults['replace'] = 'maxnonsat'
        options['replace'] = ProcessImagesPar.valid_rejection_replacements()
        dtypes['replace'] = str
        descr['replace'] = 'If all pixels are rejected, replace them using this method.  ' \
                           'Options are: {0}'.format(', '.join(options['replace']))

        defaults['lamaxiter'] = 1
        dtypes['lamaxiter'] = int
        descr['lamaxiter'] = 'Maximum number of iterations for LA cosmics routine.'

        defaults['grow'] = 1.5
        dtypes['grow'] = [int, float]
        descr['grow'] = 'Factor by which to expand regions with cosmic rays detected by the ' \
                        'LA cosmics routine.'

        defaults['rmcompact'] = True
        dtypes['rmcompact'] = bool
        descr['rmcompact'] = 'Remove compact detections in LA cosmics routine'

        defaults['sigclip'] = 4.5
        dtypes['sigclip'] = [int, float]
        descr['sigclip'] = 'Sigma level for rejection in LA cosmics routine'

        defaults['sigfrac'] = 0.3
        dtypes['sigfrac'] = [int, float]
        descr['sigfrac'] = 'Fraction for the lower clipping threshold in LA cosmics routine.'

        defaults['objlim'] = 3.0
        dtypes['objlim'] = [int, float]
        descr['objlim'] = 'Object detection limit in LA cosmics routine'

        # Instantiate the parameter set
        super(ProcessImagesPar, self).__init__(list(pars.keys()),
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
        parkeys = [ 'overscan', 'overscan_par', 'match', 'combine', 'satpix', 'sigrej', 'n_lohi',
                    'sig_lohi', 'replace', 'lamaxiter', 'grow', 'rmcompact', 'sigclip', 'sigfrac',
                    'objlim' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_overscan():
        """
        Return the valid overscan methods.
        """
        return [ 'polynomial', 'savgol', 'median' ]

    @staticmethod
    def valid_combine_methods():
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
        """
        Check the parameters are valid for the provided method.
        """

        if self.data['n_lohi'] is not None and len(self.data['n_lohi']) != 2:
            raise ValueError('n_lohi must be a list of two numbers.')
        if self.data['sig_lohi'] is not None and len(self.data['sig_lohi']) != 2:
            raise ValueError('n_lohi must be a list of two numbers.')

        if self.data['overscan'] is None:
            return
        if self.data['overscan_par'] is None:
            raise ValueError('No overscan method parameters defined!')

        # Convert param to list
        if isinstance(self.data['overscan_par'], int):
            self.data['overscan_par'] = [self.data['overscan_par']]
        
        if self.data['overscan'] == 'polynomial' and len(self.data['overscan_par']) != 3:
            raise ValueError('For polynomial overscan method, set overscan_par = order, '
                             'number of pixels, number of repeats')

        if self.data['overscan'] == 'savgol' and len(self.data['overscan_par']) != 2:
            raise ValueError('For savgol overscan method, set overscan_par = order, window size')
            
        if self.data['overscan'] == 'median' and self.data['overscan_par'] is not None:
            warnings.warn('No parameters necessary for median overscan method.  Ignoring input.')

    def to_header(self, hdr):
        """
        Write the parameters to a header object.
        """
        hdr['OSCANMET'] = (self.data['overscan'], 'Method used for overscan subtraction')
        hdr['OSCANPAR'] = (','.join([ '{0:d}'.format(p) for p in self.data['overscan_par'] ]),
                                'Overscan method parameters')
        hdr['COMBMAT'] = ('{0}'.format(self.data['match']), 'Frame combination matching')
        hdr['COMBMETH'] = (self.data['combine'], 'Method used to combine frames')
        hdr['COMBSATP'] = (self.data['satpix'], 'Saturated pixel handling when combining frames')
        hdr['COMBSIGR'] = ('{0}'.format(self.data['sigrej']),
                                'Cosmic-ray sigma rejection when combining')
        hdr['COMBNLH'] = (','.join([ '{0}'.format(n) for n in self.data['n_lohi']]),
                                'N low and high pixels rejected when combining')
        hdr['COMBSLH'] = (','.join([ '{0:.1f}'.format(s) for s in self.data['sig_lohi']]),
                                'Low and high sigma rejection when combining')
        hdr['COMBREPL'] = (self.data['replace'], 'Method used to replace pixels when combining')
        hdr['LACMAXI'] = ('{0}'.format(self.data['lamaxiter']), 'Max iterations for LA cosmic')
        hdr['LACGRW'] = ('{0:.1f}'.format(self.data['grow']), 'Growth radius for LA cosmic')
        hdr['LACRMC'] = (str(self.data['rmcompact']), 'Compact objects removed by LA cosmic')
        hdr['LACSIGC'] = ('{0:.1f}'.format(self.data['sigclip']), 'Sigma clip for LA cosmic')
        hdr['LACSIGF'] = ('{0:.1f}'.format(self.data['sigfrac']),
                            'Lower clip threshold for LA cosmic')
        hdr['LACOBJL'] = ('{0:.1f}'.format(self.data['objlim']),
                            'Object detect limit for LA cosmic')

    @classmethod
    def from_header(cls, hdr):
        """
        Instantiate the object from parameters read from a fits header.
        """
        return cls(overscan=hdr['OSCANMET'],
                   overscan_par=[int(p) for p in hdr['OSCANPAR'].split(',')],
                   match=eval(hdr['COMBMAT']),
                   combine=hdr['COMBMETH'], satpix=hdr['COMBSATP'],
                   sigrej=eval(hdr['COMBSIGR']),
                   n_lohi=[int(p) for p in hdr['COMBNLH'].split(',')],
                   sig_lohi=[float(p) for p in hdr['COMBSLH'].split(',')],
                   replace=hdr['COMBREPL'],
                   lamaxiter=int(hdr['LACMAXI']), grow=float(hdr['LACGRW']),
                   rmcompact=eval(hdr['LACRMC']), sigclip=float(hdr['LACSIGC']),
                   sigfrac=float(hdr['LACSIGF']), objlim=float(hdr['LACOBJL']))


class FlatFieldPar(ParSet):
    """
    A parameter set holding the arguments for how to perform the field
    flattening.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """
    def __init__(self, method=None, frame=None, illumflatten=None, spec_samp_fine=None, spec_samp_coarse=None,
                 spat_samp=None, tweak_slits=None, tweak_slits_thresh=None, tweak_slits_maxfrac=None):

    
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


        # ToDO there are only two methods. The bspline method and skip, so maybe we should rename the bspline method.
        defaults['method'] = 'bspline'
        options['method'] = FlatFieldPar.valid_methods()
        dtypes['method'] = str
        descr['method'] = 'Method used to flat field the data; use skip to skip flat-fielding.  ' \
                          'Options are: None, {0}'.format(', '.join(options['method']))

        # TODO: Provide a list of valid masters to use as options?
        defaults['frame'] = 'pixelflat'
        dtypes['frame'] = str
        descr['frame'] = 'Frame to use for field flattening.  Options are: "pixelflat", ' \
                         'or a specified calibration filename.'

        defaults['illumflatten'] = True
        dtypes['illumflatten'] = bool
        descr['illumflatten'] = 'Use the flat field to determine the illumination profile of each slit.'

        defaults['spec_samp_fine'] = 1.2
        dtypes['spec_samp_fine'] = [int, float]
        descr['spec_samp_fine'] = 'bspline break point spacing in units of pixels for spectral fit to flat field blaze function.'

        defaults['spec_samp_coarse'] = 50.0
        dtypes['spec_samp_coarse'] = [int, float]
        descr['spec_samp_coarse'] = 'bspline break point spacing in units of pixels for 2-d bspline-polynomial fit to ' \
                                    'flat field image residuals. This should be a large number unless you are trying to ' \
                                    'fit a sky flat with lots of narrow spectral features.'
        defaults['spat_samp'] = 5.0
        dtypes['spat_samp'] = [int, float]
        descr['spat_samp'] = 'Spatial sampling for slit illumination function. This is the width of the median ' \
                             'filter in pixels used to determine the slit illumination function, and thus sets the ' \
                             'minimum scale on which the illumination function will have features.'

        defaults['tweak_slits'] = True
        dtypes['tweak_slits'] = bool
        descr['tweak_slits'] = 'Use the illumination flat field to tweak the slit edges. ' \
                               'This will work even if illumflatten is set to False '

        defaults['tweak_slits_thresh'] = 0.93
        dtypes['tweak_slits_thresh'] = float
        descr['tweak_slits_thresh'] = 'If tweak_slits is True, this sets the illumination function threshold used to ' \
                                      'tweak the slit boundaries based on the illumination flat. ' \
                                      'It should be a number less than 1.0'

        defaults['tweak_slits_maxfrac'] = 0.10
        dtypes['tweak_slits_maxfrac'] = float
        descr['tweak_slits_maxfrac'] = 'If tweak_slit is True, this sets the maximum fractional amount (of a slits width) ' \
                                       'allowed for trimming each (i.e. left and right) slit boundary, i.e. the default is 10% ' \
                                       'which means slits would shrink or grow by at most 20% (10% on each side)'


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
        parkeys = [ 'method', 'frame', 'illumflatten', 'spec_samp_fine', 'spec_samp_coarse', 'spat_samp',
                    'tweak_slits', 'tweak_slits_thresh', 'tweak_slits_maxfrac']
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_frames():
        """
        Return the valid frame types.
        """

        # ToDO JFH So won't this fail if the user tries to provide a filename?
        return ['pixelflat'] # , 'pinhole'] disabling this for now, we don't seem to be using it. JFH

    @staticmethod
    def valid_methods():
        """
        Return the valid flat-field methods
        """
        return ['bspline', 'skip'] # [ 'PolyScan', 'bspline' ]. Same here. Not sure what PolyScan is

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        # Convert param to list
        #if isinstance(self.data['params'], int):
        #    self.data['params'] = [self.data['params']]
        
        # Check that there are the correct number of parameters
        #if self.data['method'] == 'PolyScan' and len(self.data['params']) != 3:
        #    raise ValueError('For PolyScan method, set params = order, number of '
        #                     'pixels, number of repeats')
        #if self.data['method'] == 'bspline' and len(self.data['params']) != 1:
        #    raise ValueError('For bspline method, set params = spacing (integer).')
        if self.data['frame'] in FlatFieldPar.valid_frames() or self.data['frame'] is None:
            return

        # Check the frame exists
        if not os.path.isfile(self.data['frame']):
            raise ValueError('Provided frame file name does not exist: {0}'.format(
                                self.data['frame']))

        # Check that if tweak slits is true that illumflatten is alwo true
        if self.data['tweak_slits'] and not self.data['illumflatten']:
            raise ValueError('In order to tweak slits illumflatten must be set to True')



class FlexurePar(ParSet):
    """
    A parameter set holding the arguments for how to perform the flexure
    correction.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """
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
        defaults['method'] = 'skip'
        options['method'] = FlexurePar.valid_methods()
        dtypes['method'] = str
        descr['method'] = 'Method used to correct for flexure. Use skip for no correction.  If ' \
                          'slitcen is used, the flexure correction is performed before the ' \
                          'extraction of objects (not recommended).  ' \
                          'Options are: None, {0}'.format(', '.join(options['method']))

        defaults['maxshift'] = 20
        dtypes['maxshift'] = [int, float]
        descr['maxshift'] = 'Maximum allowed flexure shift in pixels.'

        defaults['spectrum'] = os.path.join(resource_filename('pypeit', 'data/sky_spec/'), 'paranal_sky.fits')
        dtypes['spectrum'] = str
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
        return ['boxcar', 'slitcen', 'skip']

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        pass
        # TODO: This has to check both the local directory and the
        # directory in the source distribution
#        if self.data['spectrum'] is not None and not os.path.isfile(self.data['spectrum']):
#            raise ValueError('Provided archive spectrum does not exist: {0}.'.format(
#                             self.data['spectrum']))


class FluxCalibrationPar(ParSet):
    """
    A parameter set holding the arguments for how to perform the flux
    calibration.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """
    def __init__(self, balm_mask_wid=None, std_file=None, std_obj_id=None, sensfunc=None,
                 star_type=None, star_mag=None, multi_det=None, telluric=None, poly_norder=None, polycorrect=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        defaults['balm_mask_wid'] = 5.
        dtypes['balm_mask_wid'] = float
        descr['balm_mask_wid'] = 'Mask width for Balmer lines in Angstroms.'

        dtypes['std_file'] = str
        descr['std_file'] = 'Standard star file to generate sensfunc'

        dtypes['std_obj_id'] = [str, int]
        descr['std_obj_id'] = 'Specifies object in spec1d file to use as standard.' \
            ' The brightest object found is used otherwise.'

        dtypes['multi_det'] = list
        descr['multi_det'] = 'List of detector numbers to splice together for multi-detector instruments (e.g. DEIMOS)' \
                        ' They are assumed to be in order of increasing wavelength' \
                        ' And that there is *no* overlap in wavelength across detectors (might be ok if there is)'

        dtypes['sensfunc'] = str
        descr['sensfunc'] = 'FITS file that contains or will contain the sensitivity function.'

        defaults['telluric'] = False
        dtypes['telluric'] = bool
        descr['telluric'] = 'If telluric=True the code creates a sintetic standard star spectrum using the Kurucz models, ' \
            'the sens func is created setting nresln=1.5 it contains the correction for telluric lines.'

        dtypes['star_type'] = str
        descr['star_type'] = 'Spectral type of the standard star (for near-IR mainly)'

        dtypes['star_mag'] = float
        descr['star_mag'] = 'Magnitude of the standard star (for near-IR mainly)'

        defaults['poly_norder'] = 5
        dtypes['poly_norder'] = int
        descr['poly_norder'] = 'Polynomial order for sensfunc fitting'

        defaults['polycorrect'] = True
        dtypes['polycorrect'] = bool
        descr['polycorrect'] = 'Whether you want to correct the sensfunc with polynomial in the telluric and recombination line regions'

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
        parkeys = ['balm_mask_wid',  'sensfunc', 'std_file', 'std_obj_id',
                   'star_type', 'star_mag', 'multi_det', 'telluric', 'poly_norder', 'polycorrect']
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['sensfunc'] is not None and self.data['std_file'] is None and not os.path.isfile(self.data['sensfunc']):
            raise ValueError('Provided sensitivity function does not exist: {0}.'.format(
                             self.data['sensfunc']))


class ManualExtractionPar(ParSet):
    """
    DEPRECATED!!

    A parameter set holding the arguments for how to perform the
    manual extraction of a spectrum.

    A list of these objects can be included in an instance of
    :class:`ExtractObjectsPar` to perform a set of user-defined
    extractions.

    For an example of how to define a series of manual extractions in
    the pypeit input file, see :ref:`pypeit_file`.

    Args:
        frame (:obj:`str`):
            The name of the fits file for a manual extraction
        spec = List of spectral positions to hand extract
        spat = List of spatial positions to hand extract
        det = List of detectors for hand extraction. This must be a list aligned with spec and spat lists, or a single integer
             which will be used for all members of that list
        fwhm = List of FWHM for hand extraction. This must be a list aligned with spec and spat lists, or a single number which will
             be used for all members of that list'


    """
    def __init__(self, frame=None, spec=None, spat = None, det = None, fwhm = None):

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
        dtypes['frame'] = str
        descr['frame'] = 'The name of the fits file for a manual extraction'

        dtypes['spec'] = [list, float, int]
        descr['spec'] = 'List of spectral positions to hand extract '

        dtypes['spat'] = [list, float, int]
        descr['spat'] = 'List of spatial positions to hand extract '

        dtypes['det'] = [list, int]
        descr['det'] = 'List of detectors for hand extraction. This must be a list aligned with spec and spat lists, or a single integer which will be used for all members of that list.  Negative values indicated negative images.'
        dtypes['fwhm'] = [list, int,float]
        descr['fwhm'] = 'List of FWHM for hand extraction. This must be a list aligned with spec and spat lists, or a single number which will be used for all members of that list'

        # Instantiate the parameter set
        super(ManualExtractionPar, self).__init__(list(pars.keys()),
                                                  values=list(pars.values()),
                                                  dtypes=list(dtypes.values()),
                                                  descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'frame', 'spec','spat','det','fwhm']
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        pass



class ManualExtractionParOld(ParSet):
    """
    A parameter set holding the arguments for how to perform the
    manual extraction of a spectrum.

    A list of these objects can be included in an instance of
    :class:`ExtractObjectsPar` to perform a set of user-defined
    extractions.

    For an example of how to define a series of manual extractions in
    the pypeit input file, see :ref:`pypeit_file`.

    Args:
        frame (:obj:`str`):
            The name of the fits file for a manual extraction

        params (:obj:`list`):
            Parameters of the manual extraction.  For example, params =
            1,1000,500,10,10 specifies the following behavior: 1 is the
            detector number, 1000 is the spatial location that the trace
            must go through, 500 is the spectral location that the trace
            must go through, and the last two numbers (10,10) are the
            widths around the stated (spatial,spectral) location that
            should also be in the trace.'
    """
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
        dtypes['frame'] = str
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


# TODO The name of this should be changed to PypeItPar so that ReducePar can be used for what is now ScienceImagePar which
# governs the ReduceClass
class ReducePar(ParSet):
    """
    The parameter set used to hold arguments for functionality relevant
    to the overal reduction of the the data.
    
    Critically, this parameter set defines the spectrograph that was
    used to collect the data and the overall pipeline used in the
    reductions.
    
    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """
    def __init__(self, spectrograph=None, detnum=None, sortroot=None, calwin=None, scidir=None,
                 qadir=None, redux_path=None, ignore_bad_headers=None):

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
        options['spectrograph'] = ReducePar.valid_spectrographs()
        dtypes['spectrograph'] = str
        descr['spectrograph'] = 'Spectrograph that provided the data to be reduced.  ' \
                                'Options are: {0}'.format(', '.join(options['spectrograph']))

        dtypes['detnum'] = [int, list]
        descr['detnum'] = 'Restrict reduction to a list of detector indices'

        dtypes['sortroot'] = str
        descr['sortroot'] = 'A filename given to output the details of the sorted files.  If ' \
                            'None, the default is the root name of the pypeit file.  If off, ' \
                            'no output is produced.'

        # TODO: Allow this to apply to each calibration frame type
        defaults['calwin'] = 0
        dtypes['calwin']   = [int, float]
        descr['calwin'] = 'The window of time in hours to search for calibration frames for a ' \
                          'science frame'

        defaults['ignore_bad_headers'] = False
        dtypes['ignore_bad_headers'] = bool
        descr['ignore_bad_headers'] = 'Ignore bad headers (NOT recommended unless you know it is safe).'

        defaults['scidir'] = 'Science'
        dtypes['scidir'] = str
        descr['scidir'] = 'Directory relative to calling directory to write science files.'

        defaults['qadir'] = 'QA'
        dtypes['qadir'] = str
        descr['qadir'] = 'Directory relative to calling directory to write quality ' \
                         'assessment files.'

        defaults['redux_path'] = os.getcwd()
        dtypes['redux_path'] = str
        descr['redux_path'] = 'Path to folder for performing reductions.'

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
        parkeys = [ 'spectrograph', 'detnum', 'sortroot', 'calwin', 'scidir', 'qadir',
                    'redux_path', 'ignore_bad_headers']
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_spectrographs():
        # WARNING: Needs this to determine the valid spectrographs.
        # Should use pypeit.spectrographs.util.valid_spectrographs
        # instead, but it causes a circular import.  Spectrographs have
        # to be redefined here.   To fix this, spectrograph specific
        # parameter sets (like DetectorPar) and where they go needs to
        # be rethought.
        return ['gemini_gnirs','keck_deimos', 'keck_lris_blue', 'keck_lris_red', 'keck_lris_red_longonly',
                'keck_nires', 'keck_hires_red', 'keck_hires_blue', 'mmt_binospec',
                'keck_nirspec_low', 'shane_kast_blue', 'shane_kast_red', 'shane_kast_red_ret',
                'tng_dolores', 'wht_isis_blue', 'vlt_xshooter_uvb', 'vlt_xshooter_vis',
                'magellan_fire', 'magellan_mage', 'vlt_xshooter_nir', 'gemini_gmos_south_ham',
                'gemini_gmos_north_e2v', 'gemini_gmos_north_ham',
                'lbt_mods1r', 'lbt_mods1b', 'lbt_mods2r', 'lbt_mods2b', 'vlt_fors2']

    def validate(self):
        pass

    
class WavelengthSolutionPar(ParSet):
    """
    The parameter set used to hold arguments for the determination of
    wavelength solution.
    
    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """
    def __init__(self, reference=None, method=None,
                 echelle = None, ech_fix_format = None, ech_nspec_coeff = None, ech_norder_coeff = None, ech_sigrej = None,
                 lamps=None, nonlinear_counts = None,
                 sigdetect=None, fwhm=None, reid_arxiv = None, nreid_min = None, cc_thresh = None, cc_local_thresh = None,
                 nlocal_cc = None, rms_threshold=None,match_toler=None, func=None, n_first=None, n_final =None,
                 sigrej_first=None, sigrej_final=None,wv_cen=None, disp=None,numsearch=None,nfitpix=None, IDpixels=None,
                 IDwaves=None, medium=None, frame=None, nsnippet=None):
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

        # TODO JFH Does sky actually do anything?
        # TODO: Only test for 'pixel' is ever used. I.e. 'arc' or 'sky'
        # does not make a difference.
        defaults['reference'] = 'arc'
        options['reference'] = WavelengthSolutionPar.valid_reference()
        dtypes['reference'] = str
        descr['reference'] = 'Perform wavelength calibration with an arc, sky frame.  Use ' \
                             '\'pixel\' for no wavelength solution.'

        defaults['method'] = 'holy-grail'
        options['method'] = WavelengthSolutionPar.valid_methods()
        dtypes['method'] = str
        descr['method'] = 'Method to use to fit the individual arc lines. Most of these methods are now deprecated ' \
                          'as they fail most of the time without significant parameter tweaking. ' \
                          '\'holy-grail\' attempts to get a first guess at line IDs by looking for patterns in the ' \
                          'line locations. It is fully automated and works really well excpet for when it does not' \
                          '\'reidentify\' is now the preferred method, however it requires that an archive of ' \
                          'wavelength solution has been constructed for your instrument/grating combination'' \
                          ''Options are: {0}'.format(', '.join(options['method']))
#        descr['method'] = 'Method to use to fit the individual arc lines.  ' \
#                          '\'fit\' is likely more accurate, but \'simple\' uses a polynomial ' \
#                          'fit (to the log of a gaussian) and is fast and reliable.  ' \
#                          '\'arclines\' uses the arclines python package.' \
#                          'Options are: {0}'.format(', '.join(options['method']))

        # Echelle wavelength calibration stuff
        defaults['echelle'] = False
        dtypes['echelle'] = bool
        descr['echelle'] = 'Is this an echelle spectrograph? If yes an additional 2-d fit wavelength fit will be performed as a function ' \
                           'of spectral pixel and order number to improve the wavelength solution'

        defaults['ech_nspec_coeff'] = 4
        dtypes['ech_nspec_coeff'] = int
        descr['ech_nspec_coeff'] = 'For echelle spectrographs, order of the final 2d fit to the spectral dimension. ' \
                                   'You should choose this to be the n_final of the fits to the individual orders.'

        defaults['ech_norder_coeff'] = 4
        dtypes['ech_norder_coeff'] = int
        descr['ech_norder_coeff'] = 'For echelle spectrographs, order of the final 2d fit to the order dimension.'

        defaults['ech_sigrej'] = 2.0
        dtypes['ech_sigrej'] = [int,float]
        descr['ech_sigrej'] = 'For echelle spectrographs sigma clipping rejection threshold in 2d fit to spectral and order dimensions'


        # TODO: These needs to be tidied up so we can check for valid lamps. Right now I'm not checking.
        # Force lamps to be a list
        if pars['lamps'] is not None and not isinstance(pars['lamps'], list):
            pars['lamps'] = [pars['lamps']]
        options['lamps'] = None
        #options['lamps'] = WavelengthSolutionPar.valid_lamps()
        dtypes['lamps'] = list
        descr['lamps'] = 'Name of one or more ions used for the wavelength calibration.  Use ' \
                         'None for no calibration.  ' \
                         'Options are: {0}'.format(', '.join(WavelengthSolutionPar.valid_lamps()))


        # ToDo Should this be in counts or ADU? Currently the arcs are in ADU (which actually sort of makes sense here) but the
        # name of the parameter is counts. Perhaps we should just change this to nonlinear_adu or something to avoid confusion.

        # These are the parameters used for arc line detection
        # TODO: Why is this not always defined by the detectors of the
        # spectrograph?
        defaults['nonlinear_counts'] = 1e10
        dtypes['nonlinear_counts'] = float
        descr['nonlinear_counts'] = 'Arc lines above this saturation threshold are not used in wavelength solution fits because they cannot' \
                                    'be accurately centroided'

        defaults['sigdetect'] = 5.
        dtypes['sigdetect'] =  [int, float, list, numpy.ndarray]
        descr['sigdetect'] = 'Detection threshold for arc lines. This can be a single number or a list/array providing the value for each slit'

        defaults['fwhm'] = 4.
        dtypes['fwhm'] = [int, float]
        descr['fwhm'] = 'Spectral sampling of the arc lines. This is the FWHM of an arcline in *unbinned* pixels.'

        # These are the parameters used for reidentification
        defaults['reid_arxiv']=None
        dtypes['reid_arxiv'] = str
        descr['reid_arxiv'] = 'Name of the archival wavelength solution file that will be used for the wavelength ' \
                              'reidentification if the wavelength solution method = reidentify'

        defaults['nreid_min'] = 1
        dtypes['nreid_min'] = int
        descr['nreid_min'] = 'Minimum number of times that a given candidate reidentified line must be properly matched ' \
                             'with a line in the arxiv to be considered a good reidentification. If there is a lot of ' \
                             'duplication in the arxiv of the spectra in question (i.e. multislit) set this to a number ' \
                             'like 1-4. For echelle this depends on the number of solutions in the arxiv. For fixed format ' \
                             'echelle (ESI, X-SHOOTER, NIRES) set this 1. For an echelle with a tiltable grating, it will ' \
                             'depend on the number of solutions in the arxiv.'

        defaults['nsnippet'] = 2
        dtypes['nsnippet'] = int
        descr['nsnippet'] = 'Number of spectra to chop the arc spectrum into when using the full_template method'

        defaults['cc_thresh'] = 0.70
        dtypes['cc_thresh'] = [float, list, numpy.ndarray]
        descr['cc_thresh'] = 'Threshold for the *global* cross-correlation coefficient between an input spectrum and member ' \
                             'of the archive required to attempt reidentification. Spectra from the archive with a lower ' \
                             'cross-correlation are not used for reidentification. This can be a single number or a list/array providing the value for each slit'

        defaults['cc_local_thresh'] = 0.70
        dtypes['cc_local_thresh'] = float
        descr['cc_local_thresh'] = 'Threshold for the *local* cross-correlation coefficient, evaluated at each reidentified line,  ' \
                                   'between an input spectrum and the shifted and stretched archive spectrum above which a ' \
                                   'line must be to be considered a good line for reidentification. The local cross-correlation ' \
                                   'is evaluated at each candidate reidentified line (using a window of nlocal_cc), and is then ' \
                                   'used to score the the reidentified lines to arrive at the final set of good reidentifications'

        defaults['nlocal_cc'] = 11
        dtypes['nlocal_cc'] = int
        descr['nlocal_cc'] = 'Size of pixel window used for local cross-correlation computation for each arc line. If not ' \
                             'an odd number one will be added to it to make it odd.'

        defaults['ech_fix_format'] = True
        dtypes['ech_fix_format'] = bool
        descr['ech_fix_format'] = 'Is this a fixed format echelle like ESI, X-SHOOTER, or NIRES. If so reidentification ' \
                                  'will assume that each order in the data is aligned with a single order in the reid arxiv'




        # These are the parameters used for the iterative fitting of the arc lines
        defaults['rms_threshold'] = 0.15
        dtypes['rms_threshold'] = [float, list, numpy.ndarray]
        descr['rms_threshold'] = 'Minimum RMS for keeping a slit/order solution. This can be a single number or a list/array providing the value for each slit'

        defaults['match_toler'] = 2.0
        dtypes['match_toler'] = float
        descr['match_toler'] = 'Matching tolerance in pixels when searching for new lines. This is the difference ' \
                               'in pixels between the wavlength assigned to an arc line by an iteration of the wavelength ' \
                               'solution to the wavelength in the line list. This parameter is also used as the matching ' \
                               'tolerance in pixels for a line reidentification. A good line match must match within this ' \
                               'tolerance to the shifted and stretched archive spectrum, and the archive wavelength ' \
                               'solution at this match must be within match_toler dispersion elements from the line in line list.'

        defaults['func'] = 'legendre'
        dtypes['func'] = str
        descr['func'] = 'Function used for wavelength solution fits'

        defaults['n_first'] = 2
        dtypes['n_first'] = int
        descr['n_first'] = 'Order of first guess fit to the wavelength solution.'

        defaults['sigrej_first'] = 2.0
        dtypes['sigrej_first'] = float
        descr['sigrej_first'] = 'Number of sigma for rejection for the first guess to the wavelength solution.'


        defaults['n_final'] = 4
        dtypes['n_final'] = [int, float, list, numpy.ndarray]
        descr['n_final'] = 'Order of final fit to the wavelength solution. This can be a single number or a list/array providing the value for each slit'


        defaults['sigrej_final'] = 3.0
        dtypes['sigrej_final'] = float
        descr['sigrej_final'] = 'Number of sigma for rejection for the final guess to the wavelength solution.'

        # TODO: Not used
        # Backwards compatibility with basic and semi_brute algorithms
        defaults['wv_cen'] = 0.0
        dtypes['wv_cen'] = float
        descr['wv_cen'] = 'Central wavelength. Backwards compatibility with basic and semi-brute algorithms.'

        defaults['disp'] = 0.0
        dtypes['disp'] = float
        descr['disp'] = 'Dispersion. Backwards compatibility with basic and semi-brute algorithms.'


        defaults['numsearch'] = 20
        dtypes['numsearch'] = int
        descr['numsearch'] = 'Number of brightest arc lines to search for in preliminary ' \
                             'identification'

        defaults['nfitpix'] = 5
        dtypes['nfitpix'] = int
        descr['nfitpix'] = 'Number of pixels to fit when deriving the centroid of the arc ' \
                           'lines (an odd number is best)'

        dtypes['IDpixels'] = [int, float, list]
        descr['IDpixels'] = 'One or more pixels at which to manually identify a line'

        dtypes['IDwaves'] = [int, float, list]
        descr['IDwaves'] = 'Wavelengths of the manually identified lines'

        # TODO: Not used
        defaults['medium'] = 'vacuum'
        options['medium'] = WavelengthSolutionPar.valid_media()
        dtypes['medium'] = str
        descr['medium'] = 'Medium used when wavelength calibrating the data.  ' \
                          'Options are: {0}'.format(', '.join(options['medium']))

        # TODO: What should the default be?  None or 'heliocentric'?
        defaults['frame'] = 'heliocentric'
        options['frame'] = WavelengthSolutionPar.valid_reference_frames()
        dtypes['frame'] = str
        descr['frame'] = 'Frame of reference for the wavelength calibration.  ' \
                         'Options are: {0}'.format(', '.join(options['frame']))

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
        parkeys = [ 'reference', 'method',
                    'echelle', 'ech_fix_format', 'ech_nspec_coeff', 'ech_norder_coeff', 'ech_sigrej',
                    'lamps', 'nonlinear_counts', 'sigdetect', 'fwhm',
                    'reid_arxiv', 'nreid_min', 'cc_thresh', 'cc_local_thresh', 'nlocal_cc',
                    'rms_threshold', 'match_toler', 'func', 'n_first','n_final', 'sigrej_first', 'sigrej_final',
                    'wv_cen', 'disp', 'numsearch', 'nfitpix','IDpixels', 'IDwaves', 'medium', 'frame',
                    'nsnippet']
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_reference():
        """
        Return the valid wavelength solution methods.
        """
        return [ 'arc', 'sky', 'pixel' ]

    @staticmethod
    def valid_methods():
        """
        Return the valid wavelength solution methods.
        """
        return [ 'simple', 'semi-brute', 'basic','holy-grail', 'reidentify', 'full_template']

    @staticmethod
    def valid_lamps():
        """
        Return the valid lamp ions
        """
        return [ 'ArI', 'CdI', 'HgI', 'HeI', 'KrI', 'NeI', 'XeI', 'ZnI', 'ThAr' ]

    @staticmethod
    def valid_media():
        """
        Return the valid media for the wavelength calibration.
        """
        return [ 'vacuum', 'air' ]

    @staticmethod
    def valid_reference_frames():
        """
        Return the valid reference frames for the wavelength calibration
        """
        return [ 'heliocentric', 'barycentric' ]

    def validate(self):
        pass


class TraceSlitsPar(ParSet):
    """
    The parameter set used to hold arguments for tracing the slit
    positions along the dispersion axis.
    
    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """
    def __init__(self, function=None, medrep=None, number=None, trim=None,
                 maxgap=None, maxshift=None, pad=None, sigdetect=None,
                 min_slit_width = None, add_slits=None, rm_slits=None,
                 diffpolyorder=None, single=None, sobel_mode=None, pcatype=None, pcapar=None,
                 pcaextrap=None, smash_range=None, trace_npoly=None, mask_frac_thresh=None):

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
        dtypes['function'] = str
        descr['function'] = 'Function use to trace the slit center.  ' \
                            'Options are: {0}'.format(', '.join(options['function']))

        #defaults['polyorder'] = 3
        #dtypes['polyorder'] = int
        #descr['polyorder'] = 'Order of the function to use.'

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
        defaults['trim'] = (0,0)
        dtypes['trim'] = tuple
        descr['trim'] = 'How much to trim off each edge of each slit.  Each number should be 0 or positive'

        dtypes['maxgap'] = int
        descr['maxgap'] = 'Maximum number of pixels to allow for the gap between slits.  Use ' \
                          'None if the neighbouring slits are far apart or of similar ' \
                          'illumination.'

        defaults['maxshift'] = 0.15
        dtypes['maxshift'] = [int, float]
        descr['maxshift'] = 'Maximum shift in trace crude. Use a larger number for more curved slits/orders.'

        defaults['pad'] = 0
        dtypes['pad'] = int
        descr['pad'] = 'Integer number of pixels to consider beyond the slit edges.'

        defaults['sigdetect'] = 20.0
        dtypes['sigdetect'] = [int, float]
        descr['sigdetect'] = 'Sigma detection threshold for edge detection'

        defaults['mask_frac_thresh'] = 0.6
        dtypes['mask_frac_thresh'] = float
        descr['mask_frac_thresh'] = 'Minimum fraction of the slit edge that was *not* masked to use in initial PCA.'

        defaults['smash_range'] = [0., 1.]
        dtypes['smash_range'] = list
        descr['smash_range'] = 'Range of the slit in the spectral direction (in fractional units) to smash when searching for slit edges. ' \
                             'If the spectrum covers only a portion of the image, use that range.'

        defaults['trace_npoly'] = 5
        dtypes['trace_npoly'] = int
        descr['trace_npoly'] = 'Order of legendre polynomial fits to slit/order boundary traces.'

        defaults['min_slit_width'] = 6.0  # arcseconds!
        dtypes['min_slit_width'] = float
        descr['min_slit_width'] = 'If a slit spans less than this number of arcseconds over the spatial ' \
                                  'direction of the detector, it will be ignored. Use this option to prevent the ' \
                                  'of alignment (box) slits from multislit reductions, which typically cannot be reduced ' \
                                  'without a significant struggle'

        defaults['diffpolyorder'] = 2
        dtypes['diffpolyorder'] = int
        descr['diffpolyorder'] = 'Order of the 2D function used to fit the 2d solution for the ' \
                                 'spatial size of all orders.'

        # DEPRECATED
        defaults['single'] = []
        dtypes['single'] = list
        descr['single'] = 'Add a single, user-defined slit based on its location on each ' \
                          'detector.  Syntax is a list of values, 2 per detector, that define ' \
                          'the slit according to column values.  The second value (for the ' \
                          'right edge) must be greater than 0 to be applied.  LRISr example: ' \
                          'setting single = -1, -1, 7, 295 means the code will skip the ' \
                          'user-definition for the first detector but adds one for the second. ' \
                          ' None means no user-level slits defined.'

        defaults['add_slits'] = []
        dtypes['add_slits'] = [str, list]
        descr['add_slits'] = 'Add one or more user-defined slits.  This is a list of lists, with ' \
                             'each sub-list having syntax (all integers):  det:spec:spat_left:spat_right' \
                             'where det=detector, spec=spectral pixel, spat_left=spatial pixel of left slit boundary, ' \
                             ' spat_righ=spatial pixel of right slit boundary.' \
                             'For example,  2:2000:2121:2322,3:2000:1201:1500  will add a slit to detector 2 ' \
                             'passing through spec=2000 extending spatially from 2121 to 2322 and another on detector 3 ' \
                             'at spec=2000 extending from 1201 to 1500'

        defaults['rm_slits'] = []
        dtypes['rm_slits'] = [str, list]
        descr['rm_slits'] = 'Remove one or more user-specified slits.  This is a list of lists, ' \
                            'with each sub-list having syntax (all integers):  det:spec:spat ' \
                            'where det=detector, spec=spectral pixel, spat=spatial pixel.'\
                            'for example,  2:2000:2121,3:2000:1500' \
                            'the slit tracing code will remove the slits on detector 2 that contain pixel ' \
                            '(spat,spec)=(2000,2121) and detector 3 which contain (2000,2121)'

        defaults['sobel_mode'] = 'nearest'
        options['sobel_mode'] = TraceSlitsPar.valid_sobel_modes()
        dtypes['sobel_mode'] = str
        descr['sobel_mode'] = 'Mode for Sobel filtering.  Default is \'nearest\' but the ' \
                              'developers find \'constant\' works best for DEIMOS.'

        # DEPRECATED
        defaults['pcatype'] = 'pixel'
        options['pcatype'] = TraceSlitsPar.valid_pca_types()
        dtypes['pcatype'] = str
        descr['pcatype'] = 'Select to perform the PCA using the pixel position (pcatype=pixel) ' \
                           'or by spectral order (pcatype=order).  Pixel positions can be used ' \
                           'for multi-object spectroscopy where the gap between slits is ' \
                           'irregular.  Order is used for echelle spectroscopy or for slits ' \
                           'with separations that are a smooth function of the slit number.'

        # DEPRECATED
        defaults['pcapar'] = [3, 2, 1, 0]
        dtypes['pcapar'] = list
        descr['pcapar'] = 'Order of the polynomials to be used to fit the principle ' \
                          'components.  The list length must be equal to or less than ' \
                          'polyorder+1. TODO: Provide more explanation'

        defaults['pcaextrap'] = [0, 0]
        dtypes['pcaextrap'] = list
        descr['pcaextrap'] = 'The number of extra orders to predict in the negative (first ' \
                             'number) and positive (second number) direction.  Must be two ' \
                             'numbers in the list and they must be integers.'

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
        parkeys = [ 'function', 'medrep', 'number', 'trim', 'maxgap', 'maxshift',
                    'pad', 'sigdetect', 'min_slit_width', 'diffpolyorder', 'single', 'sobel_mode',
                    'pcatype', 'pcapar', 'pcaextrap', 'add_slits', 'rm_slits', 'smash_range', 'trace_npoly',
                    'mask_frac_thresh']
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
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

    @staticmethod
    def valid_pca_types():
        """
        Return the valid PCA types.
        """
        return ['pixel', 'order']

    def validate(self):
        if self.data['number'] == 0:
            raise ValueError('Number of slits must be -1 for automatic identification or '
                             'greater than 0')
        if len(self.data['pcaextrap']) != 2:
            raise ValueError('PCA extrapolation parameters must be a list with two values.')
        for e in self.data['pcaextrap']:
            if not isinstance(e, int):
                raise ValueError('PCA extrapolation values must be integers.')

class WaveTiltsPar(ParSet):
    """
    The parameter set used to hold arguments for tracing the
    monochromatic tilt along the slit.
    
    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.

    .. todo::
        Changed to reflect wavetilts.py settings.  Was `yorder`
        previously `disporder`?  If so, I think I prefer the generality
        of `disporder`...
    """
    def __init__(self, idsonly=None, tracethresh=None, sig_neigh=None, nfwhm_neigh=None, maxdev_tracefit=None, sigrej_trace=None, spat_order=None, spec_order=None,
                 func2d=None, maxdev2d=None, sigrej2d=None):


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

        #maxdev_tracefit = 1.0,
        #sigrej_trace = 3.0, max_badpix_frac = 0.20, tcrude_nave = 5,
        #npca = 1, coeff_npoly_pca = 1, sigrej_pca = 2.0,

        defaults['idsonly'] = False
        dtypes['idsonly'] = bool
        descr['idsonly'] = 'Only use the arc lines that have an identified wavelength to trace ' \
                           'tilts'

        defaults['tracethresh'] = 20.
        dtypes['tracethresh'] = [int, float, list, numpy.ndarray]
        descr['tracethresh'] = 'Significance threshold for arcs to be used in tracing wavelength tilts. ' \
                               'This can be a single number or a list/array providing the value for each slit'


        defaults['sig_neigh'] = 10.
        dtypes['sig_neigh'] = [int, float]
        descr['sig_neigh'] = 'Significance threshold for arcs to be used in line identification for the purpose of identifying neighboring lines.' \
                             'The tracethresh parameter above determines the significance threshold of lines that will be traced, but these lines' \
                             ' must be at least nfwhm_neigh fwhm away from neighboring lines. This parameter determines the significance above which' \
                             ' a line must be to be considered a possible colliding neighbor. A low value of sig_neigh will result in an overall' \
                             ' larger number of lines, which will result in more lines above tracethresh getting rejected'

        defaults['nfwhm_neigh'] = 3.0
        dtypes['nfwhm_neigh'] = [int, float]
        descr['nfwhm_neigh'] = 'Required separation between neighboring arc lines for them to be considered for tilt tracing in units of the ' \
                               'the spectral fwhm (see wavelength parset where fwhm is defined)'

        defaults['maxdev_tracefit'] = 0.2
        dtypes['maxdev_tracefit'] = [int, float]
        descr['maxdev_tracefit'] = 'Maximum absolute deviation (in units of fwhm) for the legendre polynomial fits to individual ' \
                                   'arc line tilt fits during iterative trace fitting (flux weighted, then gaussian weighted)'

        defaults['sigrej_trace'] = 3.0
        dtypes['sigrej_trace'] = [int, float]
        descr['sigrej_trace'] = 'Outlier rejection significance to determine which traced arc lines should be included in the global fit'

        defaults['spat_order'] = 3
        dtypes['spat_order'] = [int, float, list, numpy.ndarray]
        descr['spat_order'] = 'Order of the legendre polynomial to be fit to the the tilt of an arc line. This parameter determines' \
                              'both the orer of the *individual* arc line tilts, as well as the order of the spatial direction of the' \
                              '2d legendre polynomial (spatial, spectral) that is fit to obtain a global solution for the tilts across the' \
                              'slit/order. This can be a single number or a list/array providing the value for each slit'

        defaults['spec_order'] = 4
        dtypes['spec_order'] = [int, float, list, numpy.ndarray]
        descr['spec_order'] = 'Order of the spectral direction of the 2d legendre polynomial (spatial, spectral) that is ' \
                              'fit to obtain a global solution for the tilts across the slit/order. ' \
                              'This can be a single number or a list/array providing the value for each slit'


        defaults['func2d'] = 'legendre2d'
        dtypes['func2d'] = str
        descr['func2d'] = 'Type of function for 2D fit'

        defaults['maxdev2d'] = 0.25
        dtypes['maxdev2d'] = [int, float]
        descr['maxdev2d'] = 'Maximum absolute deviation (in units of fwhm) rejection threshold used to determines which pixels in global 2d fits to ' \
                            'arc line tilts are rejected because they deviate from the model by more than this value'

        defaults['sigrej2d'] = 3.0
        dtypes['sigrej2d'] = [int, float]
        descr['sigrej2d'] = 'Outlier rejection significance determining which pixels on a fit to an arc line tilt ' \
                            'are rejected by the global 2D fit'


        # Right now this is not used the fits are hard wired to be legendre for the individual fits.
        #defaults['function'] = 'legendre'
        # TODO: Allowed values?
        #dtypes['function'] = str
        #descr['function'] = 'Type of function for arc line fits'

        #defaults['yorder'] = 4
        #dtypes['yorder'] = int
        #descr['yorder'] = 'Order of the polynomial function to be used to fit the tilts ' \
        #                  'along the y direction.'


        #defaults['method'] = 'spca'
        #options['method'] = WaveTiltsPar.valid_methods()
        #dtypes['method'] = str
        #descr['method'] = 'Method used to trace the tilt of the slit along an order.  ' \
        #                  'Options are: {0}'.format(', '.join(options['method']))

        # TODO: Need to add checks that check params against method
        #defaults['params'] = [ 1, 1, 0 ]
        #dtypes['params'] = [ int, list ]
        #descr['params'] = 'Parameters to use for the provided method.  TODO: Need more explanation'

        # Instantiate the parameter set
        super(WaveTiltsPar, self).__init__(list(pars.keys()),
                                            values=list(pars.values()),
                                            defaults=list(defaults.values()),
                                            options=list(options.values()),
                                            dtypes=list(dtypes.values()),
                                            descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        parkeys = [ 'idsonly', 'tracethresh', 'sig_neigh', 'maxdev_tracefit', 'sigrej_trace','nfwhm_neigh', 'spat_order', 'spec_order', 'func2d','maxdev2d', 'sigrej2d']
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)


    def validate(self):
        pass

    #@staticmethod
    #def valid_methods():
    #    """
    #    Return the valid methods to use for tilt tracing.
    #    """
    #    return [ 'pca', 'spca', 'spline', 'interp', 'perp', 'zero' ]

#    def validate(self):
#        # Convert param to list
#        if isinstance(self.data['params'], int):
#            self.data['params'] = [self.data['params']]
#        pass

# TODO: JFH. This parameter class is now deprecated
# From artrace.trace_objects_in_slit
#       trim=2, triml=None, trimr=None, sigmin=2.0, bgreg=None
class TraceObjectsPar(ParSet):
    """
    The parameter set used to hold arguments for tracing one or more
    objects within a slit.
    
    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
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
        dtypes['function'] = str
        descr['function'] = 'Function to use to trace the object in each slit.  ' \
                            'Options are: {0}'.format(options['function'])

        defaults['order'] = 2
        dtypes['order'] = int
        descr['order'] = 'Order of the function to use to fit the object trace in each slit'

        defaults['find'] = 'standard'
        options['find'] = TraceObjectsPar.valid_detection_algorithms()
        dtypes['find'] = str
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
        dtypes['method'] = str
        descr['method'] = 'Method to use for tracing each object; only used with Echelle ' \
                          'pipeline.  Options are: {0}'.format(', '.join(options['method']))

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
        # Convert param to list
        if isinstance(self.data['params'], int):
            self.data['params'] = [self.data['params']]
        pass

# TODO JFH This parset is now deprecated.
class ExtractObjectsPar(ParSet):
    """
    The parameter set used to hold arguments for extracting object
    spectra.
    
    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
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
        dtypes['pixelmap'] = str
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
        dtypes['profile'] = str
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
        kwargs['manual'] = util.get_parset_list(cfg, 'manual', ManualExtractionPar)
        return cls(**kwargs)

    @staticmethod
    def valid_profiles():
        """
        Return the list of valid functions to use for object tracing.
        """
        return [ 'gaussian', 'gaussfunc', 'moffat', 'moffatfunc' ]

    def validate(self):
        pass

# ToDO place holder to be updated by JFH
class ScienceImagePar(ParSet):
    """
    The parameter set used to hold arguments for sky subtraction, object
    finding and extraction in the ScienceImage class

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """

    def __init__(self, bspline_spacing=None, boxcar_radius=None, trace_npoly=None,
                 global_sky_std=None, sig_thresh=None, maxnumber=None, sn_gauss=None, model_full_slit=None,
                 no_poly=None, manual=None, sky_sigrej=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k, values[k]) for k in args[1:]])  # "1:" to skip 'self'

        # Check the manual input
        if manual is not None:
            if not isinstance(manual, (ParSet, dict, list)):
                raise TypeError('Manual extraction input must be a ParSet, dictionary, or list.')
            _manual = [manual] if isinstance(manual, (ParSet, dict)) else manual
            pars['manual'] = _manual

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        defaults['bspline_spacing'] = 0.6
        dtypes['bspline_spacing'] = [int, float]
        descr['bspline_spacing'] = 'Break-point spacing for the bspline sky subtraction fits.'

        defaults['sky_sigrej'] = 3.0
        dtypes['sky_sigrej'] = float
        descr['sky_sigrej'] = 'Rejection parameter for local sky subtraction'

        defaults['boxcar_radius'] = 1.5
        dtypes['boxcar_radius'] = [int, float]
        descr['boxcar_radius'] = 'Boxcar radius in arcseconds used for boxcar extraction'

        defaults['trace_npoly'] = 5
        dtypes['trace_npoly'] = int
        descr['trace_npoly'] = 'Order of legendre polynomial fits to object traces.'

        defaults['global_sky_std'] = True
        dtypes['global_sky_std'] = bool
        descr['global_sky_std'] = 'Global sky subtraction will be performed on standard stars. This should be turned' \
                                  'off for example for near-IR reductions with narrow slits, since bright standards can' \
                                  'fill the slit causing global sky-subtraction to fail. In these situations we go ' \
                                  'straight to local sky-subtraction since it is designed to deal with such situations'

        defaults['sig_thresh'] = 10.0
        dtypes['sig_thresh'] = [int, float]
        descr['sig_thresh'] = 'Significance threshold for object finding.'

        defaults['maxnumber'] = 10
        dtypes['maxnumber'] = int
        descr['maxnumber'] = 'Maximum number of objects to extract in a science frame.  Use ' \
                             'None for no limit.'

        defaults['sn_gauss'] = 4.0
        dtypes['sn_gauss'] = [int, float]
        descr['sn_gauss'] = 'S/N threshold for performing the more sophisticated optimal extraction which performs a ' \
                            'b-spline fit to the object profile. For S/N < sn_gauss the code will simply optimal extract' \
                            'with a Gaussian with FWHM determined from the object finding.'

        defaults['model_full_slit'] = False
        dtypes['model_full_slit'] = bool
        descr['model_full_slit'] = 'If True local sky subtraction will be performed on the entire slit. If False, local sky subtraction will ' \
                            'be applied to only a restricted region around each object. This should be set to True for either multislit ' \
                            'observations using narrow slits or echelle observations with narrow slits'

        defaults['no_poly'] = False
        dtypes['no_poly'] = bool
        descr['no_poly'] = 'Turn off polynomial basis (Legendre) in global sky subtraction'

        dtypes['manual'] = list
        descr['manual'] = 'List of manual extraction parameter sets'

        # Instantiate the parameter set
        super(ScienceImagePar, self).__init__(list(pars.keys()),
                                              values=list(pars.values()),
                                              defaults=list(defaults.values()),
                                              options=list(options.values()),
                                              dtypes=list(dtypes.values()),
                                              descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        #ToDO change to updated param list
        parkeys = ['bspline_spacing', 'boxcar_radius', 'trace_npoly', 'global_sky_std', 'sig_thresh', 'maxnumber', 'sn_gauss',
                   'model_full_slit', 'no_poly', 'manual', 'sky_sigrej']
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        kwargs['manual'] = util.get_parset_list(cfg, 'manual', ManualExtractionPar)
        return cls(**kwargs)

    def validate(self):
        pass


class CalibrationsPar(ParSet):
    """
    The superset of parameters used to calibrate the science data.
    
    Note that there are specific defaults for each frame group that are
    different from the defaults of the abstracted :class:`FrameGroupPar`
    class.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """
    def __init__(self, caldir=None, reuse_masters=None, setup=None, trim=None, badpix=None,
                 biasframe=None, darkframe=None, arcframe=None, pixelflatframe=None,
                 pinholeframe=None, traceframe=None, standardframe=None, flatfield=None,
                 wavelengths=None, slits=None, tilts=None):

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
        defaults['caldir'] = 'MF'
        dtypes['caldir'] = str
        descr['caldir'] = 'Directory relative to calling directory to write master files.'

        # JFH Currently this parameter does nothing, since it is controlled by the run_pypeit argument -m
        options['reuse_masters'] = False #CalibrationsPar.allowed_master_options()
        dtypes['reuse_masters'] = bool
        descr['reuse_masters'] = 'If True PypeIt will reuse existing master frames rather than recreate them. If False, it will' \
                                 '  recreate the master frames. '
        dtypes['setup'] = str
        descr['setup'] = 'If masters=\'force\', this is the setup name to be used: e.g., ' \
                         'C_02_aa .  The detector number is ignored but the other information ' \
                         'must match the Master Frames in the master frame folder.'

        defaults['trim'] = True
        dtypes['trim'] = bool
        descr['trim'] = 'Trim the frame to isolate the data'

        defaults['badpix'] = True
        dtypes['badpix'] = bool
        descr['badpix'] = 'Make a bad pixel mask? Bias frames must be provided.'

        defaults['biasframe'] = FrameGroupPar(frametype='bias', number=5)
        dtypes['biasframe'] = [ ParSet, dict ]
        descr['biasframe'] = 'The frames and combination rules for the bias correction'

        defaults['darkframe'] = FrameGroupPar(frametype='bias', number=0)
        dtypes['darkframe'] = [ ParSet, dict ]
        descr['darkframe'] = 'The frames and combination rules for the dark-current correction'

        defaults['pixelflatframe'] = FrameGroupPar(frametype='pixelflat', number=5)
        dtypes['pixelflatframe'] = [ ParSet, dict ]
        descr['pixelflatframe'] = 'The frames and combination rules for the field flattening'

        defaults['pinholeframe'] = FrameGroupPar(frametype='pinhole', number=0)
        dtypes['pinholeframe'] = [ ParSet, dict ]
        descr['pinholeframe'] = 'The frames and combination rules for the pinholes'

        defaults['arcframe'] = FrameGroupPar(frametype='arc', number=1,
                                             process=ProcessImagesPar(sigrej=-1))
        dtypes['arcframe'] = [ ParSet, dict ]
        descr['arcframe'] = 'The frames and combination rules for the wavelength calibration'

        defaults['traceframe'] = FrameGroupPar(frametype='trace', number=3)
        dtypes['traceframe'] = [ ParSet, dict ]
        descr['traceframe'] = 'The frames and combination rules for images used for slit tracing'

        defaults['standardframe'] = FrameGroupPar(frametype='standard', number=1)
        dtypes['standardframe'] = [ ParSet, dict ]
        descr['standardframe'] = 'The frames and combination rules for the spectrophotometric ' \
                                 'standard observations'

        defaults['flatfield'] = FlatFieldPar()
        dtypes['flatfield'] = [ ParSet, dict ]
        descr['flatfield'] = 'Parameters used to set the flat-field procedure'

        defaults['wavelengths'] = WavelengthSolutionPar()
        dtypes['wavelengths'] = [ ParSet, dict ]
        descr['wavelengths'] = 'Parameters used to derive the wavelength solution'

        defaults['slits'] = TraceSlitsPar()
        dtypes['slits'] = [ ParSet, dict ]
        descr['slits'] = 'Define how the slits should be traced using the trace ?PINHOLE? frames'

        defaults['tilts'] = WaveTiltsPar()
        dtypes['tilts'] = [ ParSet, dict ]
        descr['tilts'] = 'Define how to tract the slit tilts using the trace frames'

        # Instantiate the parameter set
        super(CalibrationsPar, self).__init__(list(pars.keys()),
                                              values=list(pars.values()),
                                              defaults=list(defaults.values()),
                                              options=list(options.values()),
                                              dtypes=list(dtypes.values()),
                                              descr=list(descr.values()))
        #self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()

        # Basic keywords
        parkeys = [ 'caldir', 'reuse_masters', 'setup', 'trim', 'badpix' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None

        # Keywords that are ParSets
        pk = 'biasframe'
        kwargs[pk] = FrameGroupPar.from_dict('bias', cfg[pk]) if pk in k else None
        pk = 'darkframe'
        kwargs[pk] = FrameGroupPar.from_dict('dark', cfg[pk]) if pk in k else None
        pk = 'arcframe'
        kwargs[pk] = FrameGroupPar.from_dict('arc', cfg[pk]) if pk in k else None
        pk = 'pixelflatframe'
        kwargs[pk] = FrameGroupPar.from_dict('pixelflat', cfg[pk]) if pk in k else None
        pk = 'pinholeframe'
        kwargs[pk] = FrameGroupPar.from_dict('pinhole', cfg[pk]) if pk in k else None
        pk = 'traceframe'
        kwargs[pk] = FrameGroupPar.from_dict('trace', cfg[pk]) if pk in k else None
        pk = 'standardframe'
        kwargs[pk] = FrameGroupPar.from_dict('standard', cfg[pk]) if pk in k else None
        pk = 'flatfield'
        kwargs[pk] = FlatFieldPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'wavelengths'
        kwargs[pk] = WavelengthSolutionPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'slits'
        kwargs[pk] = TraceSlitsPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'tilts'
        kwargs[pk] = WaveTiltsPar.from_dict(cfg[pk]) if pk in k else None

        return cls(**kwargs)

    #@staticmethod
    #def allowed_master_options():
    #    """Return the allowed handling methods for the master frames."""
    #    return [ 'reuse', 'force' ]

    # TODO: Perform extensive checking that the parameters are valid for
    # the Calibrations class.  May not be necessary because validate will
    # be called for all the sub parameter sets, but this can do higher
    # level checks, if necessary.

    # JFH I'm not sure what to do about this function? Commentingo out for now.
    #def validate(self):
    #    if self.data['masters'] == 'force' \
    #            and (self.data['setup'] is None or len(self.data['setup']) == 0):
    #        raise ValueError('When forcing use of master frames, you must specify the setup to '
    #                         'be used using the \'setup\' keyword.')

#-----------------------------------------------------------------------------
# Parameters superset
class PypeItPar(ParSet):
    """
    The superset of parameters used by PypeIt.
    
    This is a single object used as a container for all the
    user-specified arguments used by PypeIt.
    
    To get the default parameters for a given spectrograph, e.g.::

        from pypeit.spectrographs.util import load_spectrograph

        spectrograph = load_spectrograph('shane_kast_blue')
        par = spectrograph.default_pypeit_par()

    If the user has a set of configuration alterations to be read from a
    pypeit file, e.g.::

        from pypeit.par.util import parse_pypeit_file
        from pypeit.spectrographs.util import load_spectrograph
        from pypeit.par import PypeItPar

        spectrograph = load_spectrograph('shane_kast_blue')
        spec_cfg_lines = spectrograph.default_pypeit_par().to_config()
        user_cfg_lines = parse_pypeit_file('myrdx.pypeit')[0]
        par = PypeItPar.from_cfg_lines(cfg_lines=spec_cfg_lines,
                                      merge_with=user_cfg_lines)

    To write the configuration of a given instance of :class:`PypeItPar`,
    use the :func:`to_config` function::
        
        par.to_config('mypypeitpar.cfg')

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`pypeitpar`.
    """
    def __init__(self, rdx=None, calibrations=None, scienceframe=None, scienceimage=None,
                 flexure=None, fluxcalib=None):

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
        defaults['rdx'] = ReducePar()
        dtypes['rdx'] = [ ParSet, dict ]
        descr['rdx'] = 'PypIt reduction rules.'

#        defaults['baseprocess'] = ProcessImagesPar()
#        dtypes['baseprocess'] = [ ParSet, dict ]
#        descr['baseprocess'] = 'Default-level parameters used when processing all images'

        defaults['calibrations'] = CalibrationsPar()
        dtypes['calibrations'] = [ ParSet, dict ]
        descr['calibrations'] = 'Parameters for the calibration algorithms'

        defaults['scienceframe'] = FrameGroupPar(frametype='science')
        dtypes['scienceframe'] = [ ParSet, dict ]
        descr['scienceframe'] = 'The frames and combination rules for the science observations'

        defaults['scienceimage'] = ScienceImagePar()
        dtypes['scienceimage'] = [ParSet, dict]
        descr['scienceimage'] = 'Parameters determining sky-subtraction, object finding, and ' \
                                'extraction'

        # Flexure is turned OFF by default
        defaults['flexure'] = FlexurePar()
        dtypes['flexure'] = [ ParSet, dict ]
        descr['flexure'] = 'Parameters used by the flexure-correction procedure.  Flexure ' \
                           'corrections are not performed by default.  To turn on, either ' \
                           'set the parameters in the \'flexure\' parameter group or set ' \
                           '\'flexure = True\' in the \'rdx\' parameter group to use the ' \
                           'default flexure-correction parameters.'

        # Flux calibration is turned OFF by default
        dtypes['fluxcalib'] = [ ParSet, dict ]
        descr['fluxcalib'] = 'Parameters used by the flux-calibration procedure.  Flux ' \
                             'calibration is not performed by default.  To turn on, either ' \
                             'set the parameters in the \'fluxcalib\' parameter group or set ' \
                             '\'fluxcalib = True\' in the \'rdx\' parameter group to use the ' \
                             'default flux-calibration parameters.'
        
        # Instantiate the parameter set
        super(PypeItPar, self).__init__(list(pars.keys()),
                                       values=list(pars.values()),
                                       defaults=list(defaults.values()),
                                       dtypes=list(dtypes.values()),
                                       descr=list(descr.values()))

        self.validate()

#    def update(self, par):
#        """
#        Update the current parameters.
#
#        Likely doesn't work because it isn't recursive ...
#        """
#        if not isinstance(par, PypeItPar):
#            raise TypeError('Parameters can only be updated using another instance of PypeItPar.')
#        self.data.update(par.data)

    @classmethod
    def from_cfg_file(cls, cfg_file=None, merge_with=None, evaluate=True):
        """
        Construct the parameter set using a configuration file.

        Note that::

            default = PypeItPar()
            nofile = PypeItPar.from_cfg_file()
            assert default.data == nofile.data, 'This should always pass.'

        Args:
            cfg_file (:obj:`str`, optional):
                The name of the configuration file that defines the
                default parameters.  This can be used to load a pypeit
                config file from a previous run that was constructed and
                output by pypeit.  This has to contain the full set of
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

        .. todo::
            Allow the user to add to the ignored strings.

        Returns:
            :class:`pypeit.par.core.PypeItPar`: The instance of the
            parameter set.
        """
        # Get the base parameters in a ConfigObj instance
        cfg = ConfigObj(PypeItPar().to_config() if cfg_file is None else cfg_file)

        # Get the list of other configuration parameters to merge it with
        _merge_with = [] if merge_with is None else \
                        ([merge_with] if isinstance(merge_with, str) else merge_with)
        merge_cfg = ConfigObj()
        for f in _merge_with:
            merge_cfg.merge(ConfigObj(f))

        # Merge with the defaults
        cfg.merge(merge_cfg)

        # Evaluate the strings if requested
        if evaluate:
            cfg = util.recursive_dict_evaluate(cfg)
        
        # Instantiate the object based on the configuration dictionary
        return cls.from_dict(cfg)

    @classmethod
    def from_cfg_lines(cls, cfg_lines=None, merge_with=None, evaluate=True):
        """
        Construct the parameter set using the list of string lines read
        from a config file.

        Note that::

            default = PypeItPar()
            nofile = PypeItPar.from_cfg_lines()
            assert default.data == nofile.data, 'This should always pass.'

        Args:
            cfg_lines (:obj:`list`, optional):
                A list of strings with lines read, or made to look like
                they are, from a configuration file.  This can be used
                to load lines from a previous run of pypeit that was
                constructed and output by pypeit.  This has to contain
                the full set of parameters, not just the subset to
                change.  For the latter, leave this as the default value
                (None) and use :arg:`merge_with` to provide a set of
                lines to merge with the defaults to construct the full
                parameter set.
            merge_with (:obj:`list`, optional):
                A list of strings with lines read, or made to look like
                they are, from a configuration file that should be
                merged with the lines provided by `cfg_lines`, or the
                default parameters.
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

        .. todo::
            Allow the user to add to the ignored strings.

        Returns:
            :class:`pypeit.par.core.PypeItPar`: The instance of the
            parameter set.
        """
        # Get the base parameters in a ConfigObj instance
        cfg = ConfigObj(PypeItPar().to_config() if cfg_lines is None else cfg_lines)
        
        # Merge in additional parameters
        if merge_with is not None:
            cfg.merge(ConfigObj(merge_with))

        # Evaluate the strings if requested
        if evaluate:
            cfg = util.recursive_dict_evaluate(cfg)
        
        # Instantiate the object based on the configuration dictionary
        return cls.from_dict(cfg)

    @classmethod
    def from_pypeit_file(cls, ifile, evaluate=True):
        """
        Construct the parameter set using a pypeit file.
        
        Args:
            ifile (str):
                Name of the pypeit file to read.  Expects to find setup
                and data blocks in the file.  See docs.
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

        .. todo::
            Allow the user to add to the ignored strings.

        Returns:
            :class:`pypeit.par.core.PypeItPar`: The instance of the
            parameter set.
        """
        # TODO: Need to include instrument-specific defaults somewhere...
        return cls.from_cfg_lines(merge_with=util.pypeit_config_lines(ifile), evaluate=evaluate)

    @classmethod
    def from_dict(cls, cfg):
        k = cfg.keys()
        kwargs = {}

        pk = 'rdx'
        kwargs[pk] = ReducePar.from_dict(cfg[pk]) if pk in k else None

        pk = 'calibrations'
        kwargs[pk] = CalibrationsPar.from_dict(cfg[pk]) if pk in k else None

        pk = 'scienceframe'
        kwargs[pk] = FrameGroupPar.from_dict('science', cfg[pk]) if pk in k else None

        pk = 'scienceimage'
        kwargs[pk] = ScienceImagePar.from_dict(cfg[pk]) if pk in k else None

        # Allow flexure to be turned on using cfg['rdx']
        pk = 'flexure'
        default = FlexurePar()
        kwargs[pk] = FlexurePar.from_dict(cfg[pk]) if pk in k else default

        # Allow flux calibration to be turned on using cfg['rdx']
        pk = 'fluxcalib'
        default = FluxCalibrationPar() \
                        if pk in cfg['rdx'].keys() and cfg['rdx']['fluxcalib'] else None
        kwargs[pk] = FluxCalibrationPar.from_dict(cfg[pk]) if pk in k else default

        if 'baseprocess' not in k:
            return cls(**kwargs)

        # Include any alterations to the basic processing of *all*
        # images
        self = cls(**kwargs)
        baseproc = ProcessImagesPar.from_dict(cfg['baseprocess'])
        self.sync_processing(baseproc)
        return self

    def sync_processing(self, proc_par):
        """
        Sync the processing of all the frame types based on the input
        ProcessImagesPar parameters.

        The parameters are merged in sequence starting from the
        parameter defaults, then including global adjustments provided
        by ``process``, and ending with the parameters that may have
        already been changed for each frame.

        This function can be used at anytime, but is most useful with
        the from_dict method where a ``baseprocess`` group can be
        supplied to change the processing parameters for all frames away
        from the defaults.

        Args:
            proc_par (:class:`ProcessImagesPar`):
                Effectively a new set of default image processing
                parameters for all frames.

        Raises:
            TypeError:
                Raised if the provided parameter set is not an instance
                of :class:`ProcessImagesPar`.
        """
        # Checks
        if not isinstance(proc_par, ProcessImagesPar):
            raise TypeError('Must provide an instance of ProcessImagesPar')
        
        # All the relevant ParSets are already ProcessImagesPar objects,
        # so we can work directly with the internal dictionaries.

        # Find the keys in the input that are different from the default
        default = ProcessImagesPar()
        base_diff = [ k for k in proc_par.keys() if default[k] != proc_par[k] ]

        # Calibration frames
        frames = [ f for f in self['calibrations'].keys() if 'frame' in f ]
        for f in frames:
            # Find the keys in self that are the same as the default
            frame_same = [ k for k in proc_par.keys() 
                            if self['calibrations'][f]['process'].data[k] == default[k] ]
            to_change = list(set(base_diff) & set(frame_same))
            for k in to_change:
                self['calibrations'][f]['process'].data[k] = proc_par[k]
            
        # Science frames
        frame_same = [ k for k in proc_par.keys() 
                            if self['scienceframe']['process'].data[k] == default[k] ]
        to_change = list(set(base_diff) & set(frame_same))
        for k in to_change:
            self['scienceframe']['process'].data[k] = proc_par[k]

    # TODO: Perform extensive checking that the parameters are valid for
    # a full run of PYPIT.  May not be necessary because validate will
    # be called for all the sub parameter sets, but this can do higher
    # level checks, if necessary.
    def validate(self):
        pass

#-----------------------------------------------------------------------------
# Instrument parameters

# TODO: This should probably get moved to spectrograph.py
class DetectorPar(ParSet):
    """
    The parameters used to define the salient properties of an
    instrument detector.

    These parameters should be *independent* of any specific use of the
    detector, and are used in the definition of the instruments served
    by PypeIt.

    To see the list of instruments served, a table with the the current
    keywords, defaults, and descriptions for the :class:`DetectorPar`
    class, and an explanation of how to define a new instrument, see
    :ref:`instruments`.
    """
    def __init__(self, dataext=None, specaxis=None, specflip=None, spatflip = None, xgap=None, ygap=None, ysize=None,
                 platescale=None, darkcurr=None, saturation=None, mincounts = None, nonlinear=None,
                 numamplifiers=None, gain=None, ronoise=None, datasec=None, oscansec=None,
                 suffix=None):

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

        # TODO: Should this be detector-specific, or camera-specific?
        defaults['specaxis'] = 0
        options['specaxis'] = [ 0, 1]
        dtypes['specaxis'] = int
        descr['specaxis'] = 'Spectra are dispersed along this axis. Allowed values are 0 ' \
                            '(first dimension for a numpy array shape) or 1 (second dimension for numpy array shape)'


        defaults['specflip'] = False
        dtypes['specflip'] = bool
        descr['specflip'] = 'If this is True then the dispersion dimension (specificed by the specaxis) will be ' \
                            'flipped so that wavelengths are always an increasing function of array index'

        defaults['spatflip'] = False
        dtypes['spatflip'] = bool
        descr['spatflip'] = 'If this is True then the spatial dimension will be ' \
                            'flipped so that blue orders for echelle spectra will appear spatially on the ' \
                            'left and wavelength will increase to the right'


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

        defaults['mincounts'] = -1000.0
        dtypes['mincounts'] = [ int, float ]
        descr['mincounts'] = 'Counts in a pixel below this value will be ignored as being unphysical'


        defaults['nonlinear'] = 0.86
        dtypes['nonlinear'] = [ int, float ]
        descr['nonlinear'] = 'Percentage of detector range which is linear (i.e. everything ' \
                             'above nonlinear*saturation will be flagged as saturated)'

        # gain, ronoise, datasec, and oscansec must be lists if there is
        # more than one amplifier
        defaults['numamplifiers'] = 1
        dtypes['numamplifiers'] = int
        descr['numamplifiers'] = 'Number of amplifiers'

        defaults['gain'] = 1.0 if pars['numamplifiers'] is None else [1.0]*pars['numamplifiers']
        dtypes['gain'] = [ int, float, list ]
        descr['gain'] = 'Inverse gain (e-/ADU). A list should be provided if a detector ' \
                        'contains more than one amplifier.'

        defaults['gain'] = 4.0 if pars['numamplifiers'] is None else [4.0]*pars['numamplifiers']
        dtypes['ronoise'] = [ int, float, list ]
        descr['ronoise'] = 'Read-out noise (e-). A list should be provided if a detector ' \
                           'contains more than one amplifier.'

        # TODO: Allow for None, such that the entire image is the data
        # section
        defaults['datasec'] = 'DATASEC' if pars['numamplifiers'] is None \
                                        else ['DATASEC']*pars['numamplifiers']
        dtypes['datasec'] = [str, list]
        descr['datasec'] = 'Either the data sections or the header keyword where the valid ' \
                           'data sections can be obtained, one per amplifier. If defined ' \
                           'explicitly should have the format of a numpy array slice'

        # TODO: Allow for None, such that there is no overscan region
        defaults['oscansec'] = 'BIASSEC' if pars['numamplifiers'] is None \
                                        else ['BIASSEC']*pars['numamplifiers']
        dtypes['oscansec'] = [str, list]
        descr['oscansec'] = 'Either the overscan section or the header keyword where the valid ' \
                            'data sections can be obtained, one per amplifier. If defined ' \
                            'explicitly should have the format of a numpy array slice'

        # TODO: Allow this to be None?
        defaults['suffix'] = ''
        dtypes['suffix'] = str
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
        parkeys = [ 'dataext', 'specaxis', 'specflip', 'spatflip','xgap', 'ygap', 'ysize', 'platescale', 'darkcurr',
                    'saturation', 'mincounts','nonlinear', 'numamplifiers', 'gain', 'ronoise', 'datasec',
                    'oscansec', 'suffix' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['numamplifiers'] > 1:
            keys = [ 'gain', 'ronoise', 'datasec', 'oscansec' ]
            dtype = [ (int, float), (int, float), str, str ]
            for i in range(len(keys)):
                if self.data[keys[i]] is None:
                    continue
                if not isinstance(self.data[keys[i]], list) \
                        or len(self.data[keys[i]]) != self.data['numamplifiers']:
                    raise ValueError('Provided {0} does not match amplifiers.'.format(keys[i]))

            for j in range(self.data['numamplifiers']):
                if self.data[keys[i]] is not None \
                        and not isinstance(self.data[keys[i]][j], dtype[i]):
                    TypeError('Incorrect type for {0}; should be {1}'.format(keys[i], dtype[i]))

# TODO: This should get moved to telescopes.py
class TelescopePar(ParSet):
    """
    The parameters used to define the salient properties of a telescope.

    These parameters should be *independent* of any specific use of the
    telescope.  They and are used by the :mod:`pypeit.telescopes` module
    to define the telescopes served by PypeIt, and kept as part of the
    :class:`pypeit.spectrographs.spectrograph.Spectrograph` definition of
    the instruments served by PypeIt.

    To see the list of instruments served, a table with the the current
    keywords, defaults, and descriptions for the :class:`TelescopePar`
    class, and an explanation of how to define a new instrument, see
    :ref:`instruments`.
    """
    def __init__(self, name=None, longitude=None, latitude=None, elevation=None, fratio=None,
                 diameter=None):

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
        defaults['name'] = 'KECK'
        options['name'] = TelescopePar.valid_telescopes()
        dtypes['name'] = str
        descr['name'] = 'Name of the telescope used to obtain the observations.  ' \
                        'Options are: {0}'.format(', '.join(options['name']))
        
        dtypes['longitude'] = [int, float]
        descr['longitude'] = 'Longitude of the telescope on Earth in degrees.'

        dtypes['latitude'] = [int, float]
        descr['latitude'] = 'Latitude of the telescope on Earth in degrees.'

        dtypes['elevation'] = [int, float]
        descr['elevation'] = 'Elevation of the telescope in m'

        dtypes['fratio'] = [int, float]
        descr['fratio'] = 'f-ratio of the telescope'

        dtypes['diameter'] = [int, float]
        descr['diameter'] = 'Diameter of the telescope in m'

        # Instantiate the parameter set
        super(TelescopePar, self).__init__(list(pars.keys()),
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
        parkeys = [ 'name', 'longitude', 'latitude', 'elevation', 'fratio', 'diameter' ]
        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_telescopes():
        """
        Return the valid telescopes.
        """
        return [ 'GEMINI-N','GEMINI-S', 'KECK', 'SHANE', 'WHT', 'APF', 'TNG', 'VLT', 'MAGELLAN', 'LBT' ]

    def validate(self):
        pass

    def platescale(self):
        r"""
        Return the platescale of the telescope in arcsec per mm.

        Calculated as

        .. math::
            p = \frac{206265}{f D},

        where :math:`f` is the f-ratio and :math:`D` is the diameter.
        If either of these is not available, the function returns
        `None`.
        """
        return None if self['fratio'] is None or self['diameter'] is None \
                else 206265/self['fratio']/self['diameter']/1e3


