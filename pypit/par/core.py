"""
Defines parameter sets used to set the behavior for core pypit
functionality.

These include:
    

    TraceSlitPar - Parameters for tracing slits
        PCAPar

    TraceTiltsPar - Parameters for tracing tilts of arc lines
        PCAPar

    TraceObjectPar - Parameter for tracing objects

    ExtractObjectPar - Parameters for extracting 1D object spectra

    TelescopePar - 

    DetectorPar

    FitsPar

    InstrumentPar

    FramePar


    - RunPar:
        - Parameters specific to a given execution of PypIt

    - ReducePar:
        - Parameters the define how PypIt should perform the reductions.
          These include general parameters and the following parameter
          subsets:
            - MastersPar: General parameters for how the master
              calibration frames should be generated and used.
            - OverscanPar: Methods used for the overscan subtraction.
            - FlatFieldPar: Methods used for field-flattening.
            - SkySubtractPar: Methods used for sky subtraction.
            - FlexurePar: Methods used for flexure correction.
            - WavelengthCalibratePar: Methods used for constructing the
              wavelength solution.
            - FluxCalibratePar: Methods used for flux calibration.



            - PixelCooPar: 

    - FrameGroupPar: Sets parameters that are used to group and combined
      frames.
        - CombineFramesPar - Parameters for combining frames

    - ArcLinesPar: Parameters used for constructing the wavelength
      solution

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
# Lower-level ParSets

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
        descr['frametype'] = 'Frame type'

        dtypes['useframe'] = str
        descr['useframe'] = 'A master calibrations file to use if it exists'

#       Should be part of the spectrograph parameters
#        dispaxis=None, 
#        defaults['dispaxis'] = 0
#        options['dispaxis'] = [0, 1]
#        dtypes['dispaxis'] = int
#        descr['dispaxis'] = 'The dispersion direction: 0 for row, 1 for column'

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

    @classmethod
    def from_dict(cls, cfg):
        return cls(frametype=cfg['frametype'], useframe=cfg['spect'], #dispaxis=cfg['dispaxis'],
                   number=cfg['number'], combine=CombinePar.from_dict(cfg['combine']))

    @staticmethod
    def valid_frame_types():
        """
        Return the list of valid frame types.
        """
        return [ 'bias', 'pixelflat', 'arc', 'pinhole', 'trace', 'standard', 'science' ]


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
        descr['match'] = 'Match frames with pixel counts that are within N-sigma of one another, '
                         'where match=N below.  If N < 0, nothing is matched.'

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
        descr['n_lohi'] = 'Number of pixels to reject at the lowest and highest ends of the '
                          'distribution; i.e., n_lohi = low, high.  Use None for no limit.'

        defaults['sig_lohi'] = [3.0, 3.0]
        dtypes['sig_lohi'] = list
        descr['sig_lohi'] = 'Sigma-clipping level at the low and high ends of the distribution; '
                            'i.e., sig_lohi = low, high.  Use None for no limit.'

        defaults['replace'] = 'maxnonsat'
        options['replace'] = CombineFramesPar.valid_rejection_replacements()
        dtypes['replace'] = str
        descr['replace'] = 'If all pixels are rejected, replace them using this method.  Option '
                           'are: {0}'.format(', '.join(options['replace']))

        # Instantiate the parameter set
        super(CombineFramesPar, self).__init__(list(pars.keys()),
                                               values=list(pars.values()),
                                               defaults=list(defaults.values()),
                                               options=list(options.values()),
                                               dtypes=list(dtypes.values()),
                                               descr=list(descr.values()))

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
        descr['method'] = 'Method to use to fit the individual arc lines.  '
                          'Options are: {0}'.format(options['method']) \
                          + '\'fit\' is likely more accurate, but \'simple\' uses a polynomial '
                            'fit (to the log of a gaussian) and is fast and reliable.  '
                            '\'arclines\' uses the arclines python package.'

        defaults['lamps'] = None
        options['lamps'] = WavelengthSolutionPar.valid_lamps()
        dtypes['lamps'] = [str, list]
        descr['lamps'] = 'Name of one or more ions used for the wavelength calibration.  Use None '
                         'for no calibration.  Options are: {0}'.format(options['lamps'])

        defaults['detection'] = 6.0
        dtypes['detection'] = [int, float]
        descr['detection'] = 'Detection threshold for arc lines (in standard deviation)'

        defaults['numsearch'] = 20
        dtypes['numsearch'] = int
        descr['numsearch'] = 'Number of brightest arc lines to search for in preliminary '
                             'identification'

        defaults['nfitpix'] = 5
        dtypes['nfitpix'] = int
        descr['nfitpix'] = 'Number of pixels to fit when deriving the centroid of the arc lines '
                           '(an odd number is best)'

        defaults['IDpixels'] = None
        dtypes['IDpixels'] = [int, list]
        descr['IDpixels'] = 'One or more pixels at which to manually identify a line'

        defaults['IDwaves'] = None
        dtypes['IDwaves'] = [int, float, list]
        descr['IDwaves'] = 'Wavelenths of the manually identified lines'

        # Instantiate the parameter set
        super(WavelenthSolutionPar, self).__init__(list(pars.keys()),
                                                   values=list(pars.values()),
                                                   defaults=list(defaults.values()),
                                                   options=list(options.values()),
                                                   dtypes=list(dtypes.values()),
                                                   descr=list(descr.values()))

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
        descr['method'] = 'Method used to fit the overscan.  '
                          'Options are: {0}'.format(options['method']) \
        
        defaults['params'] = None
        dtypes['params'] = [int, list]
        descr['params'] = 'Parameters for the overscan subtraction.  For \'polynomial\', set '
                          'params = order, number of pixels, number of repeats ; for \'savgol\', '
                          'set params = order, window size ; for \'median\', set params = None or '
                          'omit the keyword.'

        # Instantiate the parameter set
        super(OverscanPar, self).__init__(list(pars.keys()),
                                          values=list(pars.values()),
                                          defaults=list(defaults.values()),
                                          options=list(options.values()),
                                          dtypes=list(dtypes.values()),
                                          descr=list(descr.values()))

        # Check the parameters match the method requirements
        self._check_params()


    @classmethod
    def from_dict(cls, cfg):
        return cls(method=cfg['method'], params=cfg['params'])

    @staticmethod
    def valid_methods():
        """
        Return the valid overscane methods.
        """
        return [ 'polynomial', 'savgol', 'median' ]

    def _check_params(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['method'] == 'polynomial' and len(self.data['params']) != 3:
            raise ValueError('For polynomial overscan method, set params = order, number of '
                             'pixels, number of repeats')

        if self.data['method'] == 'savgol' and len(self.data['params']) != 2:
            raise ValueError('For savgol overscan method, set params = order, window size')
            
        if self.data['method'] == 'median' and self.data['params'] is not None:
            warnings.warn('No parameters necessary for median overscan method.  Ignoring input.')


class FlatFieldPar(ParSet):
    def __init__(self, frame=None, method=None, params=None, 2dpca=None):
    
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
        descr['frame'] = 'Frame to use for field flattening.  Options are: pixelflat, pinhole, '
                         'or a specified master calibration file.'

        defaults['method'] = 'bspline'
        options['method'] = FlatFieldPar.valid_methods()
        dtypes['method'] = str
        descr['method'] = 'Method used to flat field the data; use None to skip flat-fielding.  '
                          'Options are: {0}'.format(options['method'])

        defaults['params'] = 20
        dtypes['params'] = [int, list]
        descr['params'] = 'Flat-field method parameters.  For \'PolyScan\', set params = order, '
                          'numPixels, repeat ; for bspline, set params = spacing '

        # TODO:  How is 2dpca used?  Is it just another method that is
        # only used for ARMED?  Could we remove 2dpca and just add
        # another method option?
        defaults['2dpca'] = 0
        dtypes['2dpca'] = int
        descr['2dpca'] = 'Perform a simple 2D PCA on the echelle blaze fits if the value of this '
                         'argument is >1. The argument value is equal to the number of PCA '
                         'components. 0 means that no PCA will be performed.  **This is only used '
                         'with ARMED pipeline.'
    

        # Instantiate the parameter set
        super(FlatFieldPar, self).__init__(list(pars.keys()),
                                           values=list(pars.values()),
                                           defaults=list(defaults.values()),
                                           options=list(options.values()),
                                           dtypes=list(dtypes.values()),
                                           descr=list(descr.values()))

        # Check the parameters match the method requirements
        self.check_frame()


    @classmethod
    def from_dict(cls, cfg):
        return cls(frame=cfg['frame'], method=cfg['method'], params=cfg['params'],
                   2dpca=cfg['2dpca'])

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
        return [ None, 'PolyScan', 'bspline' ]


    def check_frame(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['frame'] in FlatFieldPar.valid_frames() or self.data['frame'] is None:
            return

        if not os.path.isfile(self.data['frame']):
            raise ValueError('Provided frame is not a valid frame type or existing file.')


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
        descr['method'] = 'Method used to correct for flexure. Use None for no correction.  If '
                          'slitcen is used, the flexure correction is performed before the '
                          'extraction of objects.  Options are: {0}'.format(options['method'])

        defaults['maxshift'] = 20
        dtypes['maxshift'] = [int, float]
        descr['maxshift'] = 'Maximum allowed flexure shift in pixels.'

        defaults['spectrum'] = None
        dtypes['spectrum'] = str
        descr['spectrum'] = 'Archive sky spectrum to be used for the flexure correction.'

        # Instantiate the parameter set
        super(FlexurePar, self).__init__(list(pars.keys()),
                                         values=list(pars.values()),
                                         defaults=list(defaults.values()),
                                         options=list(options.values()),
                                         dtypes=list(dtypes.values()),
                                         descr=list(descr.values()))


    @classmethod
    def from_dict(cls, cfg):
    def __init__(self, frame=None, method=None, params=None, 2dpca=None):
        return cls(frame=cfg['frame'], method=cfg['method'], params=cfg['params'],
                   2dpca=cfg['2dpca'])

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
        return [ None, 'boxcar', 'slitcen' ]


    def check_frame(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['frame'] in FlatFieldPar.valid_frames() or self.data['frame'] is None:
            return

        if not os.path.isfile(self.data['frame']):
            raise ValueError('Provided frame is not a valid frame type or existing file.')


#-----------------------------------------------------------------------------
# Top-level ParSets

class RunPar(ParSet):
    """
    Parameters specific to a given execution of PypIt.
    """
    def __init__(self, ncpus=None, calcheck=None, calwin=None, setup=None, qa=None, preponly=None,
                 stopcheck=None, useIDname=None, verbosity=None, caldir=None, scidir=None,
                 qadir=None, sort=None, overwrite=None):

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
        descr['verbosity']    = 'Level of screen output: 0 supresses all output; 1 provides '
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

        defaults['sort'] = None
        dtypes['sort'] = str
        descr['sort'] = 'File for the details of the sorted files.  If None, no output is created'

        defaults['overwrite'] = False
        dtypes['overwrite'] = bool
        descr['overwrite'] = 'Flag to overwrite any existing output files'

        # Instantiate the parameter set
        super(RunPar, self).__init__(list(pars.keys()),
                                     values=list(pars.values()),
                                     defaults=list(defaults.values()),
                                     options=list(options.values()),
                                     dtypes=list(dtypes.values()),
                                     descr=list(descr.values()),
                                     cfg_section='run', cfg_comment='Execution options')


    @classmethod
    def from_cfg_file(cls, cfg_file, merge_with=None, evaluate=True):
        # Read the configuration file
        cfg = ConfigObj(cfg_file)
        # Get the list of other files to merge it with
        _merge_with = [] if merge_with is None else \
                        ([merge_with] if isinstance(merge_with, str) else merge_with)
        for f in _merge_with:
            cfg.merge(ConfigObj(f))

        # Evaluate the strings if requested
        if evaluate:
            cfg['run'] = _recursive_dict_evaluate(cfg['run'])
        
        # Instantiate the object based on the configuration dictionary
        return cls.from_dict(cfg['run'])


    @classmethod
    def from_dict(cls, cfg):
        return cls(ncpus=cfg['ncpus'], calcheck=cfg['calcheck'], calwin=cfg['calwin'],
                   setup=cfg['setup'], qa=cfg['qa'], preponly=cfg['preponly'],
                   stopcheck=cfg['stopcheck'], useIDname=cfg['useIDname'],
                   verbosity=cfg['verbosity'], caldir=cfg['caldir'], scidir=cfg['scidir'],
                   qadir=cfg['qadir'], sortdir=cfg['sortdir'], overwrite=cfg['overwrite'])


class ReducePar(ParSet):
    """
    Parameters specific to the reduction procedures used by PypIt.
    """
    def __init__(self,

    # Spectrograph used to obtain the data to be reduced
    spectrograph = None
    # Default pipeline to use (TODO: Need to decide what None does...)
    pipeline = None
    # Restrict reduction to a single detector
    detnum = None
    # Trim the frame to isolate the data
    trim = True
    # Make a bad pixel mask? (This step requires bias frames)
    badpix = True
    # How to trace the slit center (pinhole, trace, science), you can
    # also specify a master calibrations file if it exists.
    slitcen = trace
    # How to flat field the data (trace), you can also specify a master
    # calibrations file if it exists.
    trace = trace
    # Determine the spatial slit profile
    slitprofile = False

    [[overscan]]
        # Method used to fit the overscan.  Options are: polynomial,
        # savgol, median .
        method = savgol
        # Parameters for the overscan subtraction.  For polynomial, set
        # params = order, number of pixels, number of repeats ; for
        # savgol, set params = order, window size ; for median, no
        # parameters are required (you can set params = None or omit the
        # keyword).
        params = 5, 65

    [[flatfield]]
        # Flat-field the data?
        perform = True
        # Frame to use for field flattening.  Options are: pixelflat,
        # pinhole, or a specified master calibration file.
        useframe = pixelflat
        # Method used to flat field the data: PolyScan or bspline
        method = bspline
        # Flat-field method parameters.  For PolyScan, set params =
        # order, numPixels, repeat . For bspline, set params = spacing .
        params = 20
        # Perform a simple 2D PCA on the echelle blaze fits if the value
        # of this argument is >1. The argument value is equal to the
        # number of PCA components. 0 means that no PCA will be
        # performed.  Only used with ARMED pipeline.  TODO: make this a
        # method with params?
        2dpca = 0

    [[flexure]]
        # Peform the flexure calibration
        perform = True
        # Method of flexure correction.  Options are None, 'boxcar',
        # 'slitcen'.  If 'slitcen' is used, the flexure correction is
        # performed before the extraction of objects
        method = boxcar
        # Maximum allowed flexure shift in pixels
        maxshift = 20
        # Archive sky spectrum to be used for the flexure correction
        spectrum = None

    [[calibrate]]
        # Wavelength calibrate the data?  Options are: air, vacuum, none
        wavelength = vacuum
        # Which reference frame do you want the data in (heliocentric,
        # barycentric, none)?
        refframe = heliocentric
        # Flux calibrate the data? Only used in ARMLSD pipeline.
        flux = False
        # Perform a non-linear correction.  Requires a series of
        # pixelflats of the same lamp and setup and with a variety of
        # exposure times and count rates in every pixel.
        nonlinear = False
        # YAML file with an existing calibration function
        sensfunc = None

    [[skysub]]
        # Perform the sky subtraction?
        perform = True
        # Method used for sky subtraction.  Options are: bspline.
        method = bspline
        # Method parameters.  For bspline, set params = spacing .
        params = 20

    
    ncpus=None, calcheck=None, calwin=None, setup=None, qa=None, preponly=None,
                 stopcheck=None, useIDname=None, load=None, directory=None):

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
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['spectrograph'] = 'shane_kast_blue'
        options['spectrograph']  = ReducePar.available_spectrographs()
        dtypes['spectrograph']   = str
        descr['spectrograph']    = 'Spectrograph used to obtain the data to be reduced.  ' \
                                   'Options are: {0}'.format(
                                       ', '.join(RunPar.available_spectrographs()))


    @staticmethod
    def available_spectrographs():
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

    

