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

    - TraceSlitPar: Parameters for tracing slits
        - PCAPar: PCA parameters

    - TraceTiltsPar: Parameters for tracing tilts of arc lines

    - TraceObjectsPar: Parameter for tracing objects

    - ExtractObjectPar: Parameters for extracting 1D object spectra
        - ManualExtractionPar: Parameters needed for manual extraction
          of 1D spectra

And those that define the parameters of the instrument used to collect
the data:

    - TelescopePar: Specifies the location and name of the telescope

    - DetectorPar: Specifies the properties of a given detector

    - CameraPar: Specifies the camera used for a given observation and
      contains the list of detectors for the camera.

    FitsPar


    FramePar


"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import os
import glob
import inspect

from configobj import ConfigObj

from .parset import ParSet

#-----------------------------------------------------------------------------
# Helper functions

def _pypit_root_directory():
    for p in sys.path:
        if 'PYPIT'.lower() in p.lower() \
                and len(glob.glob(os.path.join(p, 'pypit', 'pypit.py'))) == 1:
            return p
    raise OSError('Could not find PYPIT in system path.')


def _recursive_dict_evaluate(d):
    """
    will not work for lists of dicts...
    """
    for k in d.keys():
        if isinstance(d[k], dict):
           d[k] = _recursive_dict_evaluate(d[k])
        elif isinstance(d[k], list):
            try:
                d[k] = [ eval(e) for e in d[k] ]
            except:
                pass
        else:
            try:
                d[k] = eval(d[k])
            except:
                pass

    return d


#-----------------------------------------------------------------------------
# Reduction ParSets

class FrameGroupPar(ParSet):
    def __init__(self, frametype=None, useframe=None, number=None, combine=None):
        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['frametype'] = 'bias'
        options['frametype'] = FrameGroupPar.valid_frame_types()
        dtypes['frametype'] = str
        descr['frametype'] = 'Frame type.  ' \
                             'Options are: {0}'.format(', '.join(options['frametype']))

        dtypes['useframe'] = str
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
        return cls(frametype=frametype, useframe=cfg['useframe'], number=cfg['number'],
                   combine=CombineFramesPar.from_dict(cfg['combine']))

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
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['match'] = -1
        dtypes['match'] = [int, float]
        descr['match'] = 'Match frames with pixel counts that are within N-sigma of one ' \
                         'another, where match=N below.  If N < 0, nothing is matched.'

        defaults['method'] = 'mean'
        options['method'] = CombineFramesPar.valid_methods()
        dtypes['method'] = str
        descr['method'] = 'Method used to combine frames.  Options are: {0}'.format(
                                       ', '.join(options['method']))

        defaults['satpix'] = 'reject'
        options['satpix'] = CombineFramesPar.valid_saturation_handling()
        dtypes['satpix'] = str
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
        dtypes['replace'] = str
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
        return cls(match=cfg['match'], method=cfg['method'], satpix=cfg['satpix'],
                   cosmics=cfg['cosmics'], n_lohi=cfg['n_lohi'], sig_lohi=cfg['sig_lohi'],
                   replace=cfg['replace'])

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
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['method'] = 'savgol'
        options['method'] = OverscanPar.valid_methods()
        dtypes['method'] = str
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
        return cls(method=cfg['method'], params=cfg['params'])

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
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        # TODO: Provide a list of valid masters to use as options?
        defaults['frame'] = 'pixelflat'
        dtypes['frame'] = str
        descr['frame'] = 'Frame to use for field flattening.  Options are: pixelflat, pinhole, ' \
                         'or a specified master calibration file.'

        defaults['method'] = 'bspline'
        options['method'] = FlatFieldPar.valid_methods()
        dtypes['method'] = str
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
        return cls(frame=cfg['frame'], method=cfg['method'], params=cfg['params'],
                   twodpca=cfg['twodpca'])

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
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        defaults['method'] = 'boxcar'
        options['method'] = FlexurePar.valid_methods()
        dtypes['method'] = str
        descr['method'] = 'Method used to correct for flexure. Use None for no correction.  If ' \
                          'slitcen is used, the flexure correction is performed before the ' \
                          'extraction of objects.  ' \
                          'Options are: None, {0}'.format(', '.join(options['method']))

        defaults['maxshift'] = 20
        dtypes['maxshift'] = [int, float]
        descr['maxshift'] = 'Maximum allowed flexure shift in pixels.'

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
        return cls(method=cfg['method'], maxshift=cfg['maxshift'], spectrum=cfg['spectrum'])

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
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        defaults['medium'] = 'vacuum'
        options['medium'] = WavelengthCalibrationPar.valid_media()
        dtypes['medium'] = str
        descr['medium'] = 'Medium used when wavelength calibrating the data.  ' \
                          'Options are: {0}'.format(', '.join(options['medium']))

        defaults['refframe'] = 'heliocentric'
        options['refframe'] = WavelengthCalibrationPar.valid_reference_frames()
        dtypes['refframe'] = str
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
        return cls(medium=cfg['medium'], refframe=cfg['refframe'])

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
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

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

        dtypes['sensfunc'] = str
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
        return cls(flux=cfg['flux'], nonlinear=cfg['nonlinear'], sensfunc=cfg['sensfunc'])

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
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['method'] = 'bspline'
        options['method'] = SkySubtractionPar.valid_methods()
        dtypes['method'] = str
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
        return cls(method=cfg['method'], params=cfg['params'])

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
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        # TODO: Change pcatype to spacing='irregular' or
        # spacing='smooth'?
        defaults['pcatype'] = 'pixel'
        options['pcatype'] = PCAPar.valid_types()
        dtypes['pcatype'] = str
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
        return cls(pcatype=cfg['pcatype'], params=cfg['params'], extrapolate=cfg['extrapolate'])

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
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

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
        return cls(frame=cfg['frame'], params=cfg['params'])

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
        pars = dict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

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
        dtypes['caldir'] = str
        descr['caldir'] = 'Directory relative to calling directory to write master files.'
        
        defaults['scidir'] = 'Science'
        dtypes['scidir'] = str
        descr['scidir'] = 'Directory relative to calling directory to write science files.'
        
        defaults['qadir'] = 'QA'
        dtypes['qadir'] = str
        descr['qadir'] = 'Directory relative to calling directory to write qa files.'

        defaults['sortdir'] = None
        dtypes['sortdir'] = str
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
        return cls(ncpus=cfg['ncpus'], calcheck=cfg['calcheck'], calwin=cfg['calwin'],
                   setup=cfg['setup'], qa=cfg['qa'], preponly=cfg['preponly'],
                   stopcheck=cfg['stopcheck'], useIDname=cfg['useIDname'],
                   verbosity=cfg['verbosity'], caldir=cfg['caldir'], scidir=cfg['scidir'],
                   qadir=cfg['qadir'], sortdir=cfg['sortdir'], overwrite=cfg['overwrite'])

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
        pars = dict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        options['spectrograph'] = ReducePar.valid_spectrographs()
        dtypes['spectrograph'] = str
        descr['spectrograph'] = 'Spectrograph that provided the data to be reduced.  ' \
                                'Options are: {0}'.format(', '.join(options['spectrograph']))

        options['pipeline'] = ReducePar.valid_pipelines()
        dtypes['pipeline'] = str
        descr['pipeline'] = 'Pipeline options that pypit can use for reductions.  ' \
                            'Options are: {0}'.format(', '.join(options['pipeline']))

        dtypes['detnum'] = int
        descr['detnum'] = 'Restrict reduction to a single detector with this index'

        options['masters'] = ReducePar.allowed_master_options()
        dtypes['masters'] = str
        descr['masters'] = 'Treatment of master frames.  Use None to select the default ' \
                           'behavior (which is?), \'reuse\' to use any existing masters, and ' \
                           '\'force\' to __only__ use master frames.  ' \
                           'Options are: None, {0}'.format(', '.join(options['masters']))

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

        # TODO: What are the allowed center and edge trace frames?
        defaults['slit_center_frame'] = 'trace'
        dtypes['slit_center_frame'] = str
        descr['slit_center_frame'] = 'The frame that should be used to trace the slit ' \
                                     'centroid.  A master calibrations file can also be ' \
                                     'specified.'

        defaults['slit_edge_frame'] = 'trace'
        dtypes['slit_edge_frame'] = str
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
        return cls(spectrograph=cfg['spectrograph'], pipeline=cfg['pipeline'],
                   detnum=cfg['detnum'], masters=cfg['masters'], setup=cfg['setup'],
                   trim=cfg['trim'], badpix=cfg['badpix'],
                   slit_center_frame=cfg['slit_center_frame'],
                   slit_edge_frame=cfg['slit_edge_frame'],
                   overscan=OverscanPar.from_dict(cfg['overscan']),
                   flatfield=FlatFieldPar.from_dict(cfg['flatfield']),
                   flexure=FlexurePar.from_dict(cfg['flexure']),
                   wavecalib=WavelengthCalibrationPar.from_dict(cfg['wavecalib']),
                   fluxcalib=FluxCalibrationPar.from_dict(cfg['fluxcalib']),
                   skysubtract=SkySubtractionPar.from_dict(cfg['skysubtract']))

    @staticmethod
    def valid_spectrographs():
        """
        Return the list of allowed spectrographs for pypit reductions.

        .. todo::
            - Listed explicitly for now.  We can do something like what
              is currently done in the code by trolling through a file
              directory.
            - Should allow for a 'user' spectrograph
        """
        return [ 'apf_levy', 'keck_deimos', 'keck_hires', 'keck_lris_blue', 'keck_lris_red',
                 'shane_kast_blue', 'shane_kast_red', 'tng_dolores', 'wht_isis_blue' ]

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
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['method'] = 'arclines'
        options['method'] = WavelengthSolutionPar.valid_methods()
        dtypes['method'] = str
        descr['method'] = 'Method to use to fit the individual arc lines.  ' \
                          '\'fit\' is likely more accurate, but \'simple\' uses a polynomial ' \
                          'fit (to the log of a gaussian) and is fast and reliable.  ' \
                          '\'arclines\' uses the arclines python package.' \
                          'Options are: {0}'.format(', '.join(options['method']))

        options['lamps'] = WavelengthSolutionPar.valid_lamps()
        dtypes['lamps'] = [str, list]
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
        return cls(method=cfg['method'], lamps=cfg['lamps'], detection=cfg['detection'],
                   numsearch=cfg['numsearch'], nfitpix=cfg['nfitpix'], IDpixels=cfg['IDpixels'],
                   IDwaves=cfg['IDwaves'])

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
    def __init__(self, function=None, polyorder=None, medrep=None, number=None, maxgap=None,
                 pad=None, sigdetect=None, fracignore=None, diffpolyorder=None, single=None,
                 sobel_mode=None, pca=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = dict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        defaults['function'] = 'legendre'
        options['function'] = TraceSlitsPar.valid_functions()
        dtypes['function'] = str
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

        dtypes['maxgap'] = int
        descr['maxgap'] = 'Maximum number of pixels to allow for the gap between slits.  Use ' \
                          'None if the neighbouring slits are far apart or of similar ' \
                          'illumination.'

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
        dtypes['sobel_mode'] = str
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
        return cls(function=cfg['function'], polyorder=cfg['polyorder'], medrep=cfg['medrep'],
                   number=cfg['number'], maxgap=cfg['maxgap'], pad=cfg['pad'],
                   sigdetect=cfg['sigdetect'], fracignore=cfg['fracignore'],
                   diffpolyorder=cfg['diffpolyorder'], single=cfg['single'],
                   sobel_mode=cfg['sobel_mode'], pca=PCAPar.from_dict(cfg['pca']))

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
    """
    def __init__(self, idsonly=None, method=None, params=None, order=None, disporder=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = dict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        defaults['idsonly'] = False
        dtypes['idsonly'] = bool
        descr['idsonly'] = 'Only use the arc lines that have an identified wavelength to trace ' \
                           'tilts'

        defaults['method'] = 'spline'
        options['method'] = TraceTiltsPar.valid_methods()
        dtypes['method'] = str
        descr['method'] = 'Method used to trace the tilt of the slit along an order.  ' \
                          'Options are: {0}'.format(', '.join(options['method']))

        # TODO: Need to add checks that check params against method
        defaults['params'] = [ 1, 1, 0 ]
        dtypes['params'] = [ int, list ]
        descr['params'] = 'Parameters to use for the provided method.  TODO: Need more explanation'

        defaults['order'] = 1
        dtypes['order'] = int
        descr['order'] = 'Order of the polynomial function to be used for the tilt of an ' \
                         'individual arc line.  Must be 1 for eschelle data (ARMED pipeline).'

        defaults['disporder'] = 1
        dtypes['disporder'] = int
        descr['disporder'] = 'Order of the polynomial function to be used to fit the tilts ' \
                             'along the dispersion direction.  Only used by ARMED pipeline.'

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
        return cls(idsonly=cfg['idsonly'], method=cfg['method'], params=cfg['params'],
                   order=cfg['order'], disporder=cfg['disporder'])

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
        pars = dict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

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
        return cls(function=cfg['function'], order=cfg['order'], find=cfg['find'],
                   nsmooth=cfg['nsmooth'], xedge=cfg['xedge'], method=cfg['method'],
                   params=cfg['params'])

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
        pars = dict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Check the manual input
        if manual is not None:
            if not isinstance(manual, (ParSet, dict, list)):
                raise TypeError('Manual extraction input must be a ParSet, dictionary, or list.')
            _manual = [manual] if isinstance(manual, (ParSet,dict)) else manual
            pars['manual'] = _manual

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

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
        manual = []
        order = []
        for k in cfg.keys():
            if k == 'manual' and cfg[k] is None:
                continue
            if 'manual' in k:
               order += [ int(k.replace('manual','')) ]
               manual += [ ManualExtractionPar.from_dict(cfg[k]) ]
        if len(manual) > 0:
            srt = numpy.argsort(order)
            if numpy.array(order)[srt]-1 != numpy.arange(order[srt[-1]]):
                raise ValueError('Manual extraction definitions must be sequential and 1-indexed.')
            manual = (numpy.array(manual)[srt]).tolist()
        else:
            manual = None

        return cls(pixelmap=cfg['pixelmap'], pixelwidth=cfg['pixelwidth'], reuse=cfg['reuse'],
                   profile=cfg['profile'], maxnumber=cfg['maxnumber'], manual=manual)

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
class TelescopePar(ParSet):
    def __init__(self, name=None, longitude=None, latitude=None, elevation=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['name'] = 'keck'
        options['name'] = TelescopePar.valid_telescopes()
        dtypes['name'] = str
        descr['name'] = 'Name of the telescope used to obtain the observations.  ' \
                        'Options are: {0}'.format(', '.join(options['name']))
        
        defaults['longitude'] = 155.47833
        dtypes['longitude'] = [int, float]
        descr['longitude'] = 'Longitude of the telescope on Earth in degrees.'

        defaults['latitude'] = 19.82833
        dtypes['latitude'] = [int, float]
        descr['latitude'] = 'Latitude of the telescope on Earth in degrees.'

        defaults['elevation'] = 4160.0
        dtypes['elevation'] = [int, float]
        descr['elevation'] = 'Elevation of the telescope in m'

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
        return cls(name=cfg['name'], longitude=cfg['longitude'], latitude=cfg['latitude'],
                   elevation=cfg['elevation'])

    @staticmethod
    def valid_telescopes():
        """
        Return the valid overscane methods.
        """
        return [ 'keck', 'shane', 'wht', 'apf', 'tng' ]

    def validate(self):
        pass


class DetectorPar(ParSet):
    def __init__(self, frame=None, params=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['dataext'] = 0
        dtypes['dataext'] = int
        descr['dataext'] = 'Index of fits extension containing data'

        defaults['datasec'] = 'DATASEC'
        dtypes['datasec'] = str
        descr['datasec'] = 'Either the data sections or the header keyword where the valid ' \
                           'data sections can be obtained. If defined explicitly should have '
                           'the format of a numpy array slice'

        defaults['oscansec'] = 'BIASSEC'
        dtypes['oscansec'] = str
        descr['datasec'] = 'Either the overscan section or the header keyword where the valid ' \
                           'data sections can be obtained. If defined explicitly should have '
                           'the format of a numpy array slice'

        # TODO: Should this be detector- or instrument-specific?
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

        # TODO: Whate happens if this is None
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
        return cls(frame=cfg['frame'], params=cfg['params'])

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


class CameraPar(ParSet):
    """
    Parameters specific to the instrument that provided the data for
    PypIt to reduce.
    """
    def __init__(self, name=None, minexp=None, detectors=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = dict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Check the detector input
        if detectors is not None:
            if not isinstance(detectors, (ParSet, dict, list)):
                raise TypeError('Detector input must be a ParSet, dictionary, or list.')
            _detectors = [detectors] if isinstance(detectors, (ParSet,dict)) else detectors
            pars['detectors'] = _detectors

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        options['name'] = CameraPar.valid_cameras()
        dtypes['name'] = str
        descr['name'] = 'Name of the camera.  ' \
                        'Options are: {0}'.format(', '.join(options['name']))
        
        defaults['minexp'] = 1.0
        dtypes['minexp'] = [int, float]
        descr['minexp'] = 'Minimum exposure time in seconds'

        defaults['detectors'] = [ DetectorPar() ]
        dtypes['detectors'] = list
        descr['detectors'] = 'List of detectors'

        # Instantiate the parameter set
        super(CameraPar, self).__init__(list(pars.keys()),
                                        values=list(pars.values()),
                                        defaults=list(defaults.values()),
                                        options=list(options.values()),
                                        dtypes=list(dtypes.values()),
                                        descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        detectors = []
        order = []
        for k in cfg.keys():
            if k == 'det' and cfg[k] is None:
                continue
            if 'det' in k:
               order += [ int(k.replace('det','')) ]
               detectors += [ DetectorPar.from_dict(cfg[k]) ]
        if len(detectors) > 0:
            srt = numpy.argsort(order)
            if numpy.array(order)[srt]-1 != numpy.arange(order[srt[-1]]):
                raise ValueError('Manual extraction definitions must be sequential and 1-indexed.')
            detectors = (numpy.array(detectors)[srt]).tolist()
        else:
            detectors = None

        return cls(name=cfg['name'], minexp=cfg['minexp'], detectors=detectors)

    @staticmethod
    def valid_cameras():
        """
        Return the list of valid cameras.

        TODO: This should probably be telescope dependent...
        """
        return [ 'LEVY', 'DEIMOS', 'HIRES', 'LRISb', 'LRISr', 'NIRSPEC', 'KASTb', 'KASTr',
                 'DOLORES', 'ISISb' ]

    def validate(self):
        if self.data['detectors'] is None:
            raise ValueError('Must defined at least one detector!')


#-----------------------------------------------------------------------------
# Parameters superset

class PypitPar(ParSet):
    """
    The superset of all parameters used by PypIt.

    Users will likely always want to use the :func:`from_cfg_file`
    method to instantiate the parameter set, instead of using this
    instantiation function.
    """
    def __init__(self, run=None, rdx=None, bias=None, pixelflat=None, arc=None, pinhole=None,
                 trace=None, wavelengths=None, slits=None, tilts=None, objects=None, extract=None):

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
        pars = dict([(k,values[k]) for k in args[1:]])      # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['run'] = RunPar()
        dtypes['run'] = [ ParSet, dict ]
        descr['run'] = 'PypIt execution options.'

        defaults['rdx'] = ReducePar()
        dtypes['rdx'] = [ ParSet, dict ]
        descr['rdx'] = 'PypIt reduction rules.'

        defaults['bias'] = FrameGroupPar(frametype='bias')
        dtypes['bias'] = [ ParSet, dict ]
        descr['bias'] = 'The frames and combination rules for the bias correction'

        defaults['pixelflat'] = FrameGroupPar(frametype='pixelflat')
        dtypes['pixelflat'] = [ ParSet, dict ]
        descr['pixelflat'] = 'The frames and combination rules for the field flattening'

        defaults['arc'] = FrameGroupPar(frametype='arc')
        dtypes['arc'] = [ ParSet, dict ]
        descr['arc'] = 'The frames and combination rules for the wavelength calibration'

        defaults['pinhole'] = FrameGroupPar(frametype='pinhole')
        dtypes['pinhole'] = [ ParSet, dict ]
        descr['pinhole'] = 'The frames and combination rules for tracing the slit centroid'

        defaults['trace'] = FrameGroupPar(frametype='trace')
        dtypes['trace'] = [ ParSet, dict ]
        descr['trace'] = 'The frames and combination rules for tracing the slit edges'

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
    def from_cfg_file(cls, cfg_file=None, merge_with=None, evaluate=True):
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
                modifications are performed in series so the order the
                config files are listed is important.
            evaluate (:obj:`bool`, optional):
                Evaluate the values in the config object before
                assigning them in the subsequent parameter sets.  The
                parameters in the config file are *always* read as
                strings, so this should almost always be true; however,
                see the warning below.
                
        .. warning::

            When :arg:`evaluate` is true, the function runs `eval()` on
            all the entries in the ConfigObj dictionary.  This has the
            potential to go haywire if the name of a parameter
            unintentionally happens to be identical to a function that
            has been imported.  Of course, this can be useful by
            allowing one to define the function to use as a parameter,
            but it also one has to be careful the values that parameters
            should be allowed to have.

        Returns:
            `pypit.par.core.PypitPar`: The instance of the parameter set.

        """
        # Get the base parameters in a ConfigObj instance
        cfg = ConfigObj(PypitPar().to_config(None, just_lines=True)
                            if cfg_file is None else cfg_file)
        # Get the list of other files to merge it with
        _merge_with = [] if merge_with is None else \
                        ([merge_with] if isinstance(merge_with, str) else merge_with)
        for f in _merge_with:
            cfg.merge(ConfigObj(f))

        # Evaluate the strings if requested
        if evaluate:
            cfg = _recursive_dict_evaluate(cfg)
        
        # Instantiate the object based on the configuration dictionary
        return cls.from_dict(cfg)

    # TODO: What groups need to be defined?
    @classmethod
    def from_dict(cls, cfg):
        # TODO: Should allow some of the frame groups to be None,
        # meaning that, e.g., there is no cfg['pinhole']
        return cls(run=RunPar.from_dict(cfg['run']),
                   rdx=ReducePar.from_dict(cfg['rdx']),
                   bias=FrameGroupPar.from_dict('bias', cfg['bias']),
                   pixelflat=FrameGroupPar.from_dict('pixelflat', cfg['pixelflat']),
                   arc=FrameGroupPar.from_dict('arc', cfg['arc']),
                   pinhole=FrameGroupPar.from_dict('pinhole', cfg['pinhole']),
                   trace=FrameGroupPar.from_dict('trace', cfg['trace']),
                   wavelengths=WavelengthSolutionPar.from_dict(cfg['wavelengths']),
                   slits=TraceSlitsPar.from_dict(cfg['slits']),
                   tilts=TraceTiltsPar.from_dict(cfg['tilts']),
                   objects=TraceObjectsPar.from_dict(cfg['objects']),
                   extract=ExtractObjectsPar.from_dict(cfg['extract']))

    # TODO: Perform extensive checking that the parameters are valid for
    # a full run of PYPIT; this should probably just call validate() for
    # all the sub parameter sets
    def validate(self):
        pass

