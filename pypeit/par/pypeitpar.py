# encoding: utf-8
"""
Defines parameter sets used to set the behavior for core pypeit
functionality.

For more details on the full parameter hierarchy and a tabulated
description of the keywords in each parameter set, see :ref:`parameters`.

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

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
from pathlib import Path
import os
import warnings
import inspect
from IPython import embed
from collections import OrderedDict

import numpy as np

from configobj import ConfigObj

from pypeit.par.parset import ParSet
from pypeit.par import util
from pypeit.core.framematch import FrameTypeBitMask
from pypeit import msgs


def tuple_force(par):
    """
    Cast object as tuple.

    Args:
        par (object):
            Object to cast as a tuple.  Can be None; if so, the returned value
            is also None (*not* an empty tuple).  If this is already a tuple,
            this is the returned value.  If the input is a list with one tuple,
            the returned value is just the single tuple in the list (i.e, this
            does not convert the result to a tuple of one tuple).

    Returns:
        :obj:`tuple`: Casted result.
    """
    # Already has correct type
    if par is None or isinstance(par, tuple):
        return par

    # If the value is a list of one tuple, return the tuple.
    # TODO: This is a hack, and we should probably revisit how this is done.
    # The issue is that pypeit.par.util.eval_tuple always returns a list of
    # tuples, something that's required for allowing lists of detector mosaics.
    # But elements of TelluricPar are forced to be tuples.  When constructing
    # the parameters to use in a given run, the sequence of merging the
    # defaults, configuration-specific, and user-provided parameters leads to
    # converting these TelluricPar parameters into multiply nested tuples.  This
    # hook avoids that.
    if isinstance(par, list) and len(par) == 1 and isinstance(par[0], tuple):
        return par[0]

    return tuple(par)


class FrameGroupPar(ParSet):
    """
    An abstracted group of parameters that defines how specific types of
    frames should be grouped and combined.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    def __init__(self, frametype=None, useframe=None, exprng=None, process=None):
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

        # TODO: Add overscan parameters for each frame type?
        # TODO: JFH This is not documented. What are the options for useframe and what the  does it do?
        defaults['useframe'] = None
        dtypes['useframe'] = str
        descr['useframe'] = 'A calibrations file to use if it exists.'

        defaults['exprng'] = [None, None]
        dtypes['exprng'] = list
        descr['exprng'] = 'Used in identifying frames of this type.  This sets the minimum ' \
                          'and maximum allowed exposure times.  There must be two items in ' \
                          'the list.  Use None to indicate no limit; i.e., to select exposures ' \
                          'with any time greater than 30 sec, use exprng = [30, None].'

        defaults['process'] = ProcessImagesPar()
        dtypes['process'] = [ ParSet, dict ]
        descr['process'] = 'Low level parameters used for basic image processing'

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
        k = np.array([*cfg.keys()])
        parkeys = ['useframe', 'exprng']
        # TODO: cfg can contain frametype but it is ignored...
        allkeys = parkeys + ['process', 'frametype']
        badkeys = np.array([pk not in allkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for FrameGroupPar.'.format(k[badkeys]))
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
        return FrameTypeBitMask().keys()

    def validate(self):
        if len(self.data['exprng']) != 2:
            raise ValueError('exprng must be a list with two items.')


class ProcessImagesPar(ParSet):
    """
    The parameters needed to perform basic image processing.

    These parameters are primarily used by
    :class:`pypeit.processimages.ProcessImages`, the base class of many
    of the pypeit objects.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    def __init__(self, trim=None, apply_gain=None, orient=None,
                 overscan_method=None, overscan_par=None,
                 combine=None, satpix=None,
                 mask_cr=None, clip=None,
                 #cr_sigrej=None, 
                 n_lohi=None, #replace=None,
                 lamaxiter=None, grow=None,
                 comb_sigrej=None,
#                 calib_setup_and_bit=None,
                 rmcompact=None, sigclip=None, sigfrac=None, objlim=None,
                 use_biasimage=None, use_overscan=None, use_darkimage=None,
                 dark_expscale=None,
                 empirical_rn=None, shot_noise=None, noise_floor=None,
                 use_pixelflat=None, use_illumflat=None, use_specillum=None,
                 use_pattern=None, subtract_continuum=None, spat_flexure_correct=None):

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

        # Raw image fussing
        defaults['trim'] = True
        dtypes['trim'] = bool
        descr['trim'] = 'Trim the image to the detector supplied region'

        defaults['apply_gain'] = True
        dtypes['apply_gain'] = bool
        descr['apply_gain'] = 'Convert the ADUs to electrons using the detector gain'

        defaults['orient'] = True
        dtypes['orient'] = bool
        descr['orient'] = 'Orient the raw image into the PypeIt frame'

        # Bias, overscan, dark, pattern (i.e. detector "signal")
        defaults['use_biasimage'] = True
        dtypes['use_biasimage'] = bool
        descr['use_biasimage'] = 'Use a bias image.  If True, one or more must be supplied in the PypeIt file.'

        defaults['use_overscan'] = True
        dtypes['use_overscan'] = bool
        descr['use_overscan'] = 'Subtract off the overscan.  Detector *must* have one or code will crash.'

        defaults['overscan_method'] = 'savgol'
        options['overscan_method'] = ProcessImagesPar.valid_overscan_methods()
        dtypes['overscan_method'] = str
        descr['overscan_method'] = 'Method used to fit the overscan. ' \
                            'Options are: {0}'.format(', '.join(options['overscan_method']))

        defaults['overscan_par'] = [5, 65]
        dtypes['overscan_par'] = [int, list]
        descr['overscan_par'] = 'Parameters for the overscan subtraction.  For ' \
                                '\'polynomial\', set overcan_par = order, number of pixels, ' \
                                'number of repeats ; for \'savgol\', set overscan_par = ' \
                                'order, window size ; for \'median\', set overscan_par = ' \
                                'None or omit the keyword.'

        defaults['use_darkimage'] = False
        dtypes['use_darkimage'] = bool
        descr['use_darkimage'] = 'Subtract off a dark image.  If True, one or more darks must ' \
                                 'be provided.'

        defaults['dark_expscale'] = False
        dtypes['dark_expscale'] = bool
        descr['dark_expscale'] = 'If designated dark frames are used and have a different ' \
                                 'exposure time than the science frames, scale the counts by ' \
                                 'the by the ratio in the exposure times to adjust the dark ' \
                                 'counts for the difference in exposure time.  WARNING: You ' \
                                 'should always take dark frames that have the same exposure ' \
                                 'time as your science frames, so use this option with care!'

        defaults['use_pattern'] = False
        dtypes['use_pattern'] = bool
        descr['use_pattern'] = 'Subtract off a detector pattern. This pattern is assumed to be ' \
                               'sinusoidal along one direction, with a frequency that is ' \
                               'constant across the detector.'

        defaults['subtract_continuum'] = False
        dtypes['subtract_continuum'] = bool
        descr['subtract_continuum'] = 'Subtract off the continuum level from an image. This parameter should only ' \
                                      'be set to True to combine arcs with multiple different lamps. ' \
                                      'For all other cases, this parameter should probably be False.'

        defaults['empirical_rn'] = False
        dtypes['empirical_rn'] = bool
        descr['empirical_rn'] = 'If True, use the standard deviation in the overscan region to ' \
                                'measure an empirical readnoise to use in the noise model.'

        defaults['shot_noise'] = True
        dtypes['shot_noise'] = bool
        descr['shot_noise'] = 'Use the bias- and dark-subtracted image to calculate and include ' \
                              'electron count shot noise in the image processing error budget'

        defaults['noise_floor'] = 0.0
        dtypes['noise_floor'] = float
        descr['noise_floor'] = 'Impose a noise floor by adding the provided fraction of the ' \
                               'bias- and dark-subtracted electron counts to the error budget.  ' \
                               'E.g., a value of 0.01 means that the S/N of the counts in the ' \
                               'image will never be greater than 100.'

        # Flats
        defaults['use_pixelflat'] = True
        dtypes['use_pixelflat'] = bool
        descr['use_pixelflat'] = 'Use the pixel flat to make pixel-level corrections.  A ' \
                                 'pixelflat image must be provied.'

        defaults['use_illumflat'] = True
        dtypes['use_illumflat'] = bool
        descr['use_illumflat'] = 'Use the illumination flat to correct for the illumination ' \
                                 'profile of each slit.'

        defaults['use_specillum'] = False
        dtypes['use_specillum'] = bool
        descr['use_specillum'] = 'Use the relative spectral illumination profiles to correct ' \
                                 'the spectral illumination profile of each slit. This is ' \
                                 'primarily used for IFUs.  To use this, you must set ' \
                                 '``slit_illum_relative=True`` in the ``flatfield`` parameter set!'

        # Flexure
        defaults['spat_flexure_correct'] = False
        dtypes['spat_flexure_correct'] = bool
        descr['spat_flexure_correct'] = 'Correct slits, illumination flat, etc. for flexure'


        defaults['combine'] = 'mean'
        options['combine'] = ProcessImagesPar.valid_combine_methods()
        dtypes['combine'] = str
        descr['combine'] = 'Method used to combine multiple frames.  Options are: {0}'.format(
                                       ', '.join(options['combine']))

        defaults['clip'] = True
        dtypes['clip'] = bool
        descr['clip'] = 'Perform sigma clipping when combining.  Only used with combine=mean'

        defaults['comb_sigrej'] = None
        dtypes['comb_sigrej'] = float
        descr['comb_sigrej'] = 'Sigma-clipping level for when clip=True; ' \
                           'Use None for automatic limit (recommended).  '

        defaults['satpix'] = 'reject'
        options['satpix'] = ProcessImagesPar.valid_saturation_handling()
        dtypes['satpix'] = str
        descr['satpix'] = 'Handling of saturated pixels.  Options are: {0}'.format(
                                       ', '.join(options['satpix']))

        # TODO -- Make CR Parameters their own ParSet
        defaults['mask_cr'] = False
        dtypes['mask_cr'] = bool
        descr['mask_cr'] = 'Identify CRs and mask them'

#        # TODO: I don't think this is currently used; ``sigclip`` is used instead.
#        defaults['cr_sigrej'] = 20.0
#        dtypes['cr_sigrej'] = [int, float]
#        descr['cr_sigrej'] = 'Sigma level to reject cosmic rays (<= 0.0 means no CR removal)'

        defaults['n_lohi'] = [0, 0]
        dtypes['n_lohi'] = list
        descr['n_lohi'] = 'Number of pixels to reject at the lowest and highest ends of the ' \
                          'distribution; i.e., n_lohi = low, high.  Use None for no limit.'

        # TODO: I don't think this is currently used
#        defaults['replace'] = 'maxnonsat'
#        options['replace'] = ProcessImagesPar.valid_rejection_replacements()
#        dtypes['replace'] = str
#        descr['replace'] = 'If all pixels are rejected, replace them using this method.  ' \
#                           'Options are: {0}'.format(', '.join(options['replace']))

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

        # TODO: This is passed to lacosmic in
        # `pypeit.images.pypeitimage.PypeItImage.build_crmask`, *not*
        # `cr_sigrej`.
        defaults['sigclip'] = 4.5
        dtypes['sigclip'] = [int, float]
        descr['sigclip'] = 'Sigma level for rejection in LA cosmics routine'

        defaults['sigfrac'] = 0.3
        dtypes['sigfrac'] = [int, float]
        descr['sigfrac'] = 'Fraction for the lower clipping threshold in LA cosmics routine.'

        defaults['objlim'] = 3.0
        dtypes['objlim'] = [int, float]
        descr['objlim'] = 'Object detection limit in LA cosmics routine'

#        defaults['calib_setup_and_bit'] = None
#        dtypes['calib_setup_and_bit'] = str
#        descr['calib_setup_and_bit'] = 'Over-ride the calibration setup and bit, e.g. "A_7".  ' \
#                                       'Only recommended for use with quicklook.'

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
        k = np.array([*cfg.keys()])
        parkeys = ['trim', 'apply_gain', 'orient', 'use_biasimage', 'subtract_continuum', 'use_pattern',
                   'use_overscan', 'overscan_method', 'overscan_par', 'use_darkimage',
                   'dark_expscale', 'spat_flexure_correct', 'use_illumflat', 'use_specillum',
                   'empirical_rn', 'shot_noise', 'noise_floor', 'use_pixelflat', 'combine',
                   'satpix', #'calib_setup_and_bit',
                   'n_lohi', 'mask_cr',
                   'lamaxiter', 'grow', 'clip', 'comb_sigrej', 'rmcompact', 'sigclip',
                   'sigfrac', 'objlim']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for ProcessImagesPar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_overscan_methods():
        """
        Return the valid overscan methods.
        """
        return ['polynomial', 'savgol', 'median']

    @staticmethod
    def valid_combine_methods():
        """
        Return the valid methods for combining frames.
        """
        return ['median', 'mean' ]

    @staticmethod
    def valid_saturation_handling():
        """
        Return the valid approachs to handling saturated pixels.
        """
        return [ 'reject', 'force', 'nothing' ]

#    @staticmethod
#    def valid_rejection_replacements():
#        """
#        Return the valid replacement methods for rejected pixels.
#        """
#        return [ 'min', 'max', 'mean', 'median', 'weightmean', 'maxnonsat' ]

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """

        if self.data['n_lohi'] is not None and len(self.data['n_lohi']) != 2:
            raise ValueError('n_lohi must be a list of two numbers.')

        if not self.data['use_overscan']:
            return
        if self.data['overscan_par'] is None:
            raise ValueError('No overscan method parameters defined!')

        # Convert param to list
        if isinstance(self.data['overscan_par'], int):
            self.data['overscan_par'] = [self.data['overscan_par']]

        if self.data['overscan_method'] == 'polynomial' and len(self.data['overscan_par']) != 3:
            raise ValueError('For polynomial overscan method, set overscan_par = order, '
                             'number of pixels, number of repeats')

        if self.data['overscan_method'] == 'savgol' and len(self.data['overscan_par']) != 2:
            raise ValueError('For savgol overscan method, set overscan_par = order, window size')

        if self.data['overscan_method'] == 'median' and self.data['overscan_par'] is not None:
            warnings.warn('No parameters necessary for median overscan method.  Ignoring input.')

        if not self.data['use_pixelflat'] \
                and (self.data['use_illumflat'] or self.data['use_specillum']):
            raise ValueError('To apply a slit-illumination or spectral flat-field correction, '
                             'you must also apply the pixel-flat correction.')

    # TODO: Are these out of date or is this a purposeful subselection of the
    # full parameter set?
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
        hdr['COMBSRJ'] = (self.data['comb_sigrej'], 'Sigma rejection when combining')
#        hdr['COMBREPL'] = (self.data['replace'], 'Method used to replace pixels when combining')
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
                   n_lohi=[int(p) for p in hdr['COMBNLH'].split(',')],
                   comb_sigrej=float(hdr['COMBSRJ']),
#                   replace=hdr['COMBREPL'],
#                   cr_sigrej=eval(hdr['LASIGR']),
                   lamaxiter=int(hdr['LACMAXI']), grow=float(hdr['LACGRW']),
                   rmcompact=eval(hdr['LACRMC']), sigclip=float(hdr['LACSIGC']),
                   sigfrac=float(hdr['LACSIGF']), objlim=float(hdr['LACOBJL']))


class FlatFieldPar(ParSet):
    """
    A parameter set holding the arguments for how to perform the field
    flattening.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    def __init__(self, method=None, pixelflat_file=None, spec_samp_fine=None,
                 spec_samp_coarse=None, spat_samp=None, tweak_slits=None, tweak_slits_thresh=None,
                 tweak_slits_maxfrac=None, rej_sticky=None, slit_trim=None, slit_illum_pad=None,
                 illum_iter=None, illum_rej=None, twod_fit_npoly=None, saturated_slits=None,
                 slit_illum_relative=None, slit_illum_ref_idx=None, slit_illum_smooth_npix=None,
                 pixelflat_min_wave=None, pixelflat_max_wave=None, slit_illum_finecorr=None,
                 fit_2d_det_response=None):

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

        # Pixel flat parameter keys
        defaults['pixelflat_file'] = None
        dtypes['pixelflat_file'] = str
        descr['pixelflat_file'] = 'Filename of the image to use for pixel-level field flattening'

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

        defaults['pixelflat_min_wave'] = None
        dtypes['pixelflat_min_wave'] = [int, float]
        descr['pixelflat_min_wave'] = 'All values of the normalized pixel flat are set to 1 for wavelengths below this value.'

        defaults['pixelflat_max_wave'] = None
        dtypes['pixelflat_max_wave'] = [int, float]
        descr['pixelflat_max_wave'] = 'All values of the normalized pixel flat are set to 1 for wavelengths above this value.'


        # Slits
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


        defaults['rej_sticky'] = False
        dtypes['rej_sticky'] = bool
        descr['rej_sticky'] = 'Propagate the rejected pixels through the stages of the ' \
                              'flat-field fitting (i.e, from the spectral fit, to the spatial ' \
                              'fit, and finally to the 2D residual fit).  If False, pixels ' \
                              'rejected in each stage are included in each subsequent stage.'

        defaults['slit_trim'] = 3.
        dtypes['slit_trim'] = [int, float, tuple]
        descr['slit_trim'] = 'The number of pixels to trim each side of the slit when ' \
                             'selecting pixels to use for fitting the spectral response ' \
                             'function.  Single values are used for both slit edges; a ' \
                             'two-tuple can be used to trim the left and right sides differently.'

        defaults['slit_illum_pad'] = 5.
        dtypes['slit_illum_pad'] = [int, float]
        descr['slit_illum_pad'] = 'The number of pixels to pad the slit edges when constructing ' \
                                  'the slit-illumination profile. Single value applied to both ' \
                                  'edges.'

        defaults['slit_illum_finecorr'] = True
        dtypes['slit_illum_finecorr'] = bool
        descr['slit_illum_finecorr'] = 'If True, a fine correction to the spatial illumination profile ' \
                                       'will be performed. The fine correction is a low order 2D polynomial ' \
                                       'fit to account for a gradual change to the spatial illumination ' \
                                       'profile as a function of wavelength.'

        defaults['slit_illum_relative'] = False
        dtypes['slit_illum_relative'] = bool
        descr['slit_illum_relative'] = 'Generate an image of the relative spectral illumination ' \
                                       'for a multi-slit setup.  If you set ``use_slitillum = ' \
                                       'True`` for any of the frames that use the flatfield ' \
                                       'model, this *must* be set to True. Currently, this is ' \
                                       'only used for IFU reductions.'

        defaults['illum_iter'] = 0
        dtypes['illum_iter'] = int
        descr['illum_iter'] = 'The number of rejection iterations to perform when constructing ' \
                              'the slit-illumination profile.  No rejection iterations are ' \
                              'performed if 0.  WARNING: Functionality still being tested.'

        defaults['illum_rej'] = 5.
        dtypes['illum_rej'] = [int, float]
        descr['illum_rej'] = 'The sigma threshold used in the rejection iterations used to ' \
                             'refine the slit-illumination profile.  Rejection iterations are ' \
                             'only performed if ``illum_iter > 0``.'

        dtypes['twod_fit_npoly'] = int
        descr['twod_fit_npoly'] = 'Order of polynomial used in the 2D bspline-polynomial fit to ' \
                                  'flat-field image residuals. The code determines the order of ' \
                                  'these polynomials to each slit automatically depending on ' \
                                  'the slit width, which is why the default is None. Alter ' \
                                  'this paramter at your own risk!'

        defaults['saturated_slits'] = 'crash'
        options['saturated_slits'] = FlatFieldPar.valid_saturated_slits_methods()
        dtypes['saturated_slits'] = str
        descr['saturated_slits'] = 'Behavior when a slit is encountered with a large fraction ' \
                                   'of saturated pixels in the flat-field.  The options are: ' \
                                   '\'crash\' - Raise an error and halt the data reduction; ' \
                                   '\'mask\' - Mask the slit, meaning no science data will be ' \
                                   'extracted from the slit; \'continue\' - ignore the ' \
                                   'flat-field correction, but continue with the reduction.'

        defaults['slit_illum_ref_idx'] = 0
        dtypes['slit_illum_ref_idx'] = int
        descr['slit_illum_ref_idx'] = 'The index of a reference slit (0-indexed) used for estimating ' \
                                      'the relative spectral sensitivity (or the relative blaze). This ' \
                                      'parameter is only used if ``slit_illum_relative = True``.'

        defaults['slit_illum_smooth_npix'] = 10
        dtypes['slit_illum_smooth_npix'] = int
        descr['slit_illum_smooth_npix'] = 'The number of pixels used to determine smoothly varying ' \
                                          'relative weights is given by ``nspec/slit_illum_smooth_npix``, ' \
                                          'where nspec is the number of spectral pixels.'

        defaults['fit_2d_det_response'] = False
        dtypes['fit_2d_det_response'] = bool
        descr['fit_2d_det_response'] = 'Set this variable to True if you want to compute and ' \
                                       'account for the detector response in the flatfield image. ' \
                                       'Note that ``detector response`` refers to pixel sensitivity ' \
                                       'variations that primarily depend on (x,y) detector coordinates. ' \
                                       'In most cases, the default 2D bspline is sufficient to account ' \
                                       'for detector response (i.e. set this parameter to False). Note ' \
                                       'that this correction will _only_ be performed for the spectrographs ' \
                                       'that have a dedicated response correction implemented. Currently,' \
                                       'this correction is only implemented for Keck+KCWI.'

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
        k = np.array([*cfg.keys()])
        parkeys = ['method', 'pixelflat_file', 'spec_samp_fine', 'spec_samp_coarse',
                   'spat_samp', 'pixelflat_min_wave', 'pixelflat_max_wave',
                   'tweak_slits', 'tweak_slits_thresh', 'tweak_slits_maxfrac',
                   'rej_sticky', 'slit_trim', 'slit_illum_pad', 'slit_illum_relative',
                   'illum_iter', 'illum_rej', 'twod_fit_npoly', 'saturated_slits',
                   'slit_illum_ref_idx', 'slit_illum_smooth_npix', 'slit_illum_finecorr', 'fit_2d_det_response']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for FlatFieldPar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_methods():
        """
        Return the valid flat-field methods
        """
        return ['bspline', 'skip'] # [ 'PolyScan', 'bspline' ]. Same here. Not sure what PolyScan is

    @staticmethod
    def valid_saturated_slits_methods():
        """
        Return the valid options for dealing with saturated slits.
        """
        return ['crash', 'mask', 'continue']

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
        if self.data['pixelflat_file'] is None:
            return

        # Check the frame exists
        if not os.path.isfile(self.data['pixelflat_file']):
            raise ValueError('Provided frame file name does not exist: {0}'.format(
                                self.data['pixelflat_file']))

        # Check that if tweak slits is true that illumflatten is alwo true
        # TODO -- We don't need this set, do we??   See the desc of tweak_slits above
        #if self.data['tweak_slits'] and not self.data['illumflatten']:
        #    raise ValueError('In order to tweak slits illumflatten must be set to True')


class FlexurePar(ParSet):
    """
    A parameter set holding the arguments for how to perform flexure
    corrections, both spatial and spectral

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    def __init__(self, spec_method=None, spec_maxshift=None, spectrum=None,
                 multi_min_SN=None, excessive_shift=None):

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
        defaults['spec_method'] = 'skip'
        options['spec_method'] = FlexurePar.valid_methods()
        dtypes['spec_method'] = str
        descr['spec_method'] = 'Method used to correct for flexure. Use skip for no correction.  If ' \
                          'slitcen is used, the flexure correction is performed before the ' \
                          'extraction of objects (not recommended).  ' \
                          'Options are: None, {0}'.format(', '.join(options['spec_method']))

        defaults['spec_maxshift'] = 20
        dtypes['spec_maxshift'] = int
        descr['spec_maxshift'] = 'Maximum allowed spectral flexure shift in pixels.'

        defaults['spectrum'] = 'paranal_sky.fits'
        dtypes['spectrum'] = str
        descr['spectrum'] = 'Archive sky spectrum to be used for the flexure correction.'

        # The following are all for MultiDet flexure
        defaults['multi_min_SN'] = 1
        dtypes['multi_min_SN'] = [int, float]
        descr['multi_min_SN'] = 'Minimum S/N for analyzing sky spectrum for flexure'

        defaults['excessive_shift'] = 'use_median'
        options['excessive_shift'] = FlexurePar.valid_excessive_shift_methods()
        dtypes['excessive_shift'] = str
        descr['excessive_shift'] = 'Behavior when the measured spectral flexure shift is ' \
                                   'larger than ``spec_maxshift``.  The options are: ' \
                                   '\'crash\' - Raise an error and halt the data reduction; ' \
                                   '\'set_to_zero\' - Set the flexure shift to zero and continue ' \
                                   'with the reduction; \'continue\' - Use the large ' \
                                   'flexure value whilst issuing a warning; and \'use_median\' - ' \
                                   'Use the median flexure shift among all the objects in the same slit ' \
                                   '(if more than one object is detected) or among all ' \
                                   'the other slits; if not available, the flexure correction will not be applied.'

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

        k = np.array([*cfg.keys()])
        parkeys = ['spec_method', 'spec_maxshift', 'spectrum',
                   'multi_min_SN', 'excessive_shift']
#                   'spat_frametypes']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for FlexurePar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_methods():
        """
        Return the valid flat-field methods
        """
        return ['boxcar', 'slitcen', 'skip']

    @staticmethod
    def valid_excessive_shift_methods():
        """
        Return the valid options for dealing with excessive flexure shift.
        """
        return ['crash', 'set_to_zero', 'continue', 'use_median']

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


class AlignPar(ParSet):
    """
    The parameter set used to hold arguments for tracing the
    alignments in an align frame.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """

    def __init__(self, locations=None, trace_npoly=None, trim_edge=None, snr_thresh=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k, values[k]) for k in args[1:]])  # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        defaults['locations'] = [0.0, 0.5, 1.0]
        dtypes['locations'] = [list, np.ndarray]
        descr['locations'] = 'Locations of the bars, in a list, specified as a fraction of the slit width'

        defaults['trace_npoly'] = 4
        dtypes['trace_npoly'] = int
        descr['trace_npoly'] = 'Order of the polynomial to use when fitting the trace of a single bar'

        defaults['trim_edge'] = [0, 0]
        dtypes['trim_edge'] = list
        descr['trim_edge'] = 'Trim the slit by this number of pixels left/right before finding alignment bars'

        defaults['snr_thresh'] = 1.0  # This must be low, because the routine will find the
        dtypes['snr_thresh'] = [int, float]
        descr['snr_thresh'] = 'S/N ratio threshold for finding an alignment trace. This should be a low ' \
                              'number to ensure that the algorithm finds all bars. The algorithm will ' \
                              'then only use the N most significant detections, where N is the number ' \
                              'of elements specified in the "locations" keyword argument'

        # Instantiate the parameter set
        super(AlignPar, self).__init__(list(pars.keys()),
                                            values=list(pars.values()),
                                            defaults=list(defaults.values()),
                                            options=list(options.values()),
                                            dtypes=list(dtypes.values()),
                                            descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])
        parkeys = ['locations', 'trace_npoly', 'trim_edge', 'snr_thresh']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for WaveTiltsPar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        pass


class Coadd1DPar(ParSet):
    """
    A parameter set holding the arguments for how to perform 1D coadds

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    def __init__(self, ex_value=None, flux_value=None, nmaskedge=None,
                 sn_smooth_npix=None, wave_method=None, dv=None, wave_grid_min=None, wave_grid_max=None, spec_samp_fact=None, ref_percentile=None, maxiter_scale=None,
                 sigrej_scale=None, scale_method=None, sn_min_medscale=None, sn_min_polyscale=None, maxiter_reject=None,
                 lower=None, upper=None, maxrej=None, sn_clip=None, nbest=None, sensfuncfile=None, coaddfile=None,
                 mag_type=None, filter=None, filter_mag=None, filter_mask=None, chk_version=None):

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

        # Extraction to use
        defaults['ex_value'] = 'OPT'
        options['ex_value'] = Coadd1DPar.valid_ex()
        dtypes['ex_value'] = str
        descr['ex_value'] = "The extraction to coadd, i.e. optimal or boxcar. Must be either 'OPT' or 'BOX'"

        # Fluxed?
        defaults['flux_value'] = True
        dtypes['flux_value'] = bool
        descr['flux_value'] = 'If True (default), the code will coadd the fluxed spectra (i.e. the FLAM) in the ' \
                              'spec1d files. If False, it will coadd the counts.'

        # Mask edge pixels?
        defaults['nmaskedge'] = 2
        dtypes['nmaskedge'] = int
        descr['nmaskedge'] = 'Number of edge pixels to mask. This should be removed/fixed.'

        # Offsets
        defaults['sn_smooth_npix'] = None
        dtypes['sn_smooth_npix'] = [int, float]
        descr['sn_smooth_npix'] = 'Number of pixels to median filter by when computing S/N used to decide how to scale ' \
                                  'and weight spectra. If set to None (default), the code will determine the effective ' \
                                  'number of good pixels per spectrum in the stack that is being co-added and use 10% of ' \
                                  'this neff.'


        # Offsets
        defaults['wave_method'] = 'linear'
        dtypes['wave_method'] = str
        descr['wave_method'] = "Method used to construct wavelength grid for coadding spectra. The routine that creates " \
                               "the wavelength is :func:`~pypeit.core.wavecal.wvutils.get_wave_grid`. The options are:" \
                               " "\
                               "'iref' -- Use the first wavelength array.  " \
                               "'velocity' -- Grid is uniform in velocity.  " \
                               "'log10' -- Grid is uniform in log10(wave). This is the same as velocity.  " \
                               "'linear' -- Grid is uniform in lambda.  " \
                               "'concatenate' -- Meld the input wavelength arrays"

        defaults['dv'] = None
        dtypes['dv'] = [int, float]
        descr['dv'] = "Dispersion in units of km/s in case you want to specify it in the get_wave_grid  (for the 'velocity' option), " \
                    "otherwise a median value is computed from the data."

        defaults['wave_grid_min'] = None
        dtypes['wave_grid_min'] = [int, float]
        descr['wave_grid_min'] = "Used in case you want to specify the minimum wavelength in your wavelength grid, default=None computes from data"

        defaults['wave_grid_max'] = None
        dtypes['wave_grid_max'] = [int, float]
        descr['wave_grid_max'] = "Used in case you want to specify the maximum wavelength in your wavelength grid, default=None computes from data"

        defaults['spec_samp_fact'] = 1.0
        dtypes['spec_samp_fact'] = float
        descr['spec_samp_fact'] = "Make the wavelength grid  sampling finer (spec_samp_fact < 1.0) or coarser " \
                                  "(spec_samp_fact > 1.0) by this sampling factor. This basically multiples the 'native' " \
                                  "spectral pixels by spec_samp_fact, i.e. units spec_samp_fact are pixels."

        defaults['ref_percentile'] = 70.0
        dtypes['ref_percentile'] = [int, float]
        descr['ref_percentile'] = 'Percentile used for selecting the minimum SNR cut from a reference spectrum used to ' \
                                  'robustly determine the median ratio between spectra. This parameter is used by ' \
                                  'coadd1d.robust_median_ratio as part of the automatic rescaling procedure. Pixels ' \
                                  'above this percentile cut are deemed the "good" pixels and are used to compute the ' \
                                  'ratio of two spectra.  This must be a number between 0 and 100.'

        defaults['maxiter_scale'] = 5
        dtypes['maxiter_scale'] = int
        descr['maxiter_scale'] = 'Maximum number of iterations performed for rescaling spectra.'

        defaults['sigrej_scale'] = 3.0
        dtypes['sigrej_scale'] = [int, float]
        descr['sigrej_scale'] = 'Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec.'

        defaults['scale_method'] = 'auto'
        dtypes['scale_method'] = str
        descr['scale_method'] = "Method used to rescale the spectra prior to coadding. The options are:" \
                                " "\
                                "'auto' -- Determine the scaling method automatically based on the S/N ratio which works well.  "\
                                "'poly' -- Polynomial rescaling.  " \
                                "'median' -- Median rescaling  " \
                                "'none' -- Do not rescale.  " \
                                "'hand' -- Pass in hand scaling factors. This option is not well tested."

        defaults['sn_min_medscale'] = 0.5
        dtypes['sn_min_medscale'] = [int, float]
        descr['sn_min_medscale'] = "For scale method set to ``auto``, this sets the minimum SNR for which median scaling is attempted."

        defaults['sn_min_polyscale'] = 2.0
        dtypes['sn_min_polyscale'] = [int, float]
        descr['sn_min_polyscale'] = "For scale method set to ``auto``, this sets the minimum SNR for which polynomial scaling is attempted."


        defaults['maxiter_reject'] = 5
        dtypes['maxiter_reject'] = int
        descr['maxiter_reject'] = 'Maximum number of iterations for stacking and rejection. The code stops iterating ' \
                                  'either when the output mask does not change betweeen successive iterations or when ' \
                                  'maxiter_reject is reached.'
        defaults['lower'] = 3.0
        dtypes['lower'] = [int, float]
        descr['lower'] = 'Lower rejection threshold used for rejecting pixels when combining spectra in units of sigma.'

        defaults['upper'] = 3.0
        dtypes['upper'] = [int, float]
        descr['upper'] = 'Upper rejection threshold used for rejecting pixels when combining spectra in units of sigma.'

        defaults['maxrej'] = None
        dtypes['maxrej'] = int
        descr['maxrej'] = 'Coadding performs iterative rejection by comparing each exposure to a preliminary stack of ' \
                          'all the exposures. If this parameter is set then it will not reject more than maxrej pixels ' \
                          'per iteration of this rejection. The default is None, which means no maximum on rejected pixels.'

        defaults['sn_clip'] = 30.0
        dtypes['sn_clip'] = [int, float]
        descr['sn_clip'] = 'Errors are capped during rejection so that the S/N is never greater than sn_clip. This ' \
                           'prevents overly aggressive rejection in high S/N ratio spectrum which neverthless differ ' \
                           'at a level greater than the formal S/N due to systematics.'

        defaults['nbest'] = None
        dtypes['nbest'] = int
        descr['nbest'] = 'Number of orders to use for estimating the per exposure weights. Default is None, ' \
                         'which will just use one fourth of the total number of orders. This is only used for Echelle.'

        # For scaling to an input filter magnitude
        defaults['filter'] = 'none'
        dtypes['filter'] = str
        descr['filter'] = 'Filter for scaling.  See flux_calib.load_fitler_file() for naming.  Ignore if none'

        defaults['mag_type'] = 'AB'
        dtypes['mag_type'] = str
        descr['mag_type'] = 'Magnitude type.  AB is the only option currently allowed'

        defaults['filter_mag'] = None
        dtypes['filter_mag'] = float
        descr['filter_mag'] = 'Magnitude of the source in the given filter'

        defaults['filter_mask'] = None
        dtypes['filter_mask'] = [str, list]
        descr['filter_mask'] = 'List of wavelength regions to mask when doing the scaling (`i.e.`, occasional junk pixels). '\
                               'Colon and comma separateed, e.g.   5552:5559,6010:6030'


        # JFH These last two are actually arguments and not parameters that are only here because there is no other easy
        # way to parse .coadd1d files except with parsets. I would like to separate arguments from parameters.
        defaults['sensfuncfile'] = None
        dtypes['sensfuncfile'] = str
        descr['sensfuncfile'] = 'File containing sensitivity function which is a requirement for echelle coadds. ' \
                            'This is only used for Echelle.'

        defaults['coaddfile'] = None
        dtypes['coaddfile'] = str
        descr['coaddfile'] = 'Output filename'


        # Version checking
        defaults['chk_version'] = True
        dtypes['chk_version'] = bool
        descr['chk_version'] = 'If True enforce strict PypeIt version checking to ensure that spec1d*.fits files were created' \
                               'with the current version of PypeIt'

        # Instantiate the parameter set
        super(Coadd1DPar, self).__init__(list(pars.keys()),
                                         values=list(pars.values()),
                                         defaults=list(defaults.values()),
                                         dtypes=list(dtypes.values()),
                                         descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])
        parkeys = ['ex_value', 'flux_value', 'nmaskedge', 'sn_smooth_npix', 'wave_method', 'dv', 'wave_grid_min', 'wave_grid_max',
                   'spec_samp_fact', 'ref_percentile', 'maxiter_scale', 'sigrej_scale', 'scale_method',
                   'sn_min_medscale', 'sn_min_polyscale', 'maxiter_reject', 'lower', 'upper',
                   'maxrej', 'sn_clip', 'nbest', 'sensfuncfile', 'coaddfile', 'chk_version',
                   'filter', 'mag_type', 'filter_mag', 'filter_mask']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for Coadd1DPar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        pass

    @staticmethod
    def valid_ex():
        """
        Return the valid flat-field methods
        """
        return ['BOX', 'OPT']


class Coadd2DPar(ParSet):
    """
    A parameter set holding the arguments for how to perform 2D coadds

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    def __init__(self, only_slits=None, offsets=None, spat_toler=None, weights=None, user_obj=None,
                 use_slits4wvgrid=None, manual=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Offsets
        defaults['only_slits'] = None
        dtypes['only_slits'] = [int, list]
        descr['only_slits'] = 'Slit ID, or list of slit IDs that the user want to restrict the coadd to. ' \
                              'I.e., only this/these slit/s will be coadded.'

        defaults['offsets'] = 'auto'
        dtypes['offsets'] = [str, list]
        descr['offsets'] = 'Offsets for the images being combined (spat pixels). Options are: ' \
                           '``maskdef_offsets``, ``header``, ``auto``, and a list of offsets. ' \
                           'Use ``maskdef_offsets`` to use the offsets computed during the slitmask ' \
                           'design matching (currently available for DEIMOS and MOSFIRE only). If equal ' \
                           'to ``header``, the dither offsets recorded in the header, when available, will be used. ' \
                           'If ``auto`` is chosen, PypeIt will try to compute the offsets using a reference object ' \
                           'with the highest S/N, or an object selected by the user (see ``user_obj``). ' \
                           'If a list of offsets is provided, PypeIt will use it.'

        defaults['spat_toler'] = 5
        dtypes['spat_toler'] = int
        descr['spat_toler'] = 'This parameter provides the desired tolerance in spatial pixel used ' \
                              'to identify slits in different exposures'

        # Offsets
        defaults['use_slits4wvgrid'] = False
        dtypes['use_slits4wvgrid'] = bool
        descr['use_slits4wvgrid'] = 'If True, use the slits to set the trace down the center'

        # Weights
        defaults['weights'] = 'auto'
        dtypes['weights'] = [str, list]
        descr['weights'] = 'Mode for the weights used to coadd images. Options are: ' \
                           '``auto``, ``uniform``, or a list of weights. If ``auto`` is used, ' \
                           'PypeIt will try to compute the weights using a reference object ' \
                           'with the highest S/N, or an object selected by the user (see ``user_obj``), ' \
                           'if ``uniform`` is used, uniform weights will be applied. If a list of weights ' \
                           'is provided, PypeIt will use it.'

        # object to use for weights and offsets
        defaults['user_obj'] = None
        dtypes['user_obj'] = [int, list]
        descr['user_obj'] = 'Object that the user wants to use to compute the weights and/or the ' \
                            'offsets for coadding images. For slit spectroscopy, provide the ' \
                            '``SLITID`` and the ``OBJID``, separated by comma, of the selected object. ' \
                            'For echelle spectroscopy, provide the ``ECH_OBJID`` of the selected object. ' \
                            'See :doc:`out_spec1D` for more info about ``SLITID``, ``OBJID`` and ``ECH_OBJID``. ' \
                            'If this parameter is not ``None``, it will be used to compute the offsets ' \
                            'only if ``offsets = auto``, and it will used to compute ' \
                            'the weights only if ``weights = auto``.'

        # manual extraction
        defaults['manual'] = None
        dtypes['manual'] = str
        descr['manual'] = 'Manual extraction parameters. det:spat:spec:fwhm:boxcar_radius. ' \
                          'Multiple manual extractions are semi-colon separated, ' \
                          'and spat,spec are in the pseudo-image generated by COADD2D.' \
                              'boxcar_radius is optional and in pixels (not arcsec!).'

        # Instantiate the parameter set
        super(Coadd2DPar, self).__init__(list(pars.keys()),
                                                 values=list(pars.values()),
                                                 defaults=list(defaults.values()),
                                                 dtypes=list(dtypes.values()),
                                                 descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])
        parkeys = ['only_slits', 'offsets', 'spat_toler', 'weights', 'user_obj', 'use_slits4wvgrid', 'manual']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for Coadd2DPar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        pass


class CubePar(ParSet):
    """
    The parameter set used to hold arguments for functionality relevant
    to cube generation (primarily for IFU data).

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """

    def __init__(self, slit_spec=None, relative_weights=None, combine=None, output_filename=None,
                 standard_cube=None, reference_image=None, save_whitelight=None, method=None,
                 ra_min=None, ra_max=None, dec_min=None, dec_max=None, wave_min=None, wave_max=None,
                 spatial_delta=None, wave_delta=None, astrometric=None, grating_corr=None, scale_corr=None,
                 skysub_frame=None, spec_subpixel=None, spat_subpixel=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k, values[k]) for k in args[1:]])  # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        # Cube Parameters
        defaults['slit_spec'] = True
        dtypes['slit_spec'] = [bool]
        descr['slit_spec'] = 'If the data use slits in one spatial direction, set this to True. ' \
                             'If the data uses fibres for all spaxels, set this to False.'

        defaults['relative_weights'] = False
        dtypes['relative_weights'] = [bool]
        descr['relative_weights'] = 'If set to True, the combined frames will use a relative weighting scheme. ' \
                                    'This only works well if there is a common continuum source in the field of ' \
                                    'view of all input observations, and is generally only required if high ' \
                                    'relative precision is desired.'

        defaults['combine'] = False
        dtypes['combine'] = [bool]
        descr['combine'] = 'If set to True, the input frames will be combined. Otherwise, a separate ' \
                           'datacube will be generated for each input spec2d file, and will be saved as ' \
                           'a spec3d file.'

        defaults['output_filename'] = ""
        dtypes['output_filename'] = str
        descr['output_filename'] = 'If combining multiple frames, this string sets the output filename of ' \
                                   'the combined datacube. If combine=False, the output filenames will be ' \
                                   'prefixed with "spec3d_*"'

        defaults['standard_cube'] = None
        dtypes['standard_cube'] = str
        descr['standard_cube'] = 'Filename of a standard star datacube. This cube will be used to correct ' \
                                 'the relative scales of the slits, and to flux calibrate the science ' \
                                 'datacube.'

        defaults['reference_image'] = None
        dtypes['reference_image'] = str
        descr['reference_image'] = 'White light image of a previously combined datacube. The white light ' \
                                   'image will be used as a reference when calculating the offsets of the ' \
                                   'input spec2d files. Ideally, the reference image should have the same ' \
                                   'shape as the data to be combined (i.e. set the ra_min, ra_max etc. params ' \
                                   'so they are identical to the reference image).'

        defaults['save_whitelight'] = False
        dtypes['save_whitelight'] = bool
        descr['save_whitelight'] = 'Save a white light image of the combined datacube. The output filename ' \
                                   'will be given by the "output_filename" variable with a suffix "_whitelight". ' \
                                   'Note that the white light image collapses the flux along the wavelength axis, ' \
                                   'so some spaxels in the 2D white light image may have different wavelength ' \
                                   'ranges. If combine=False, the individual spec3d files will have a suffix ' \
                                   '"_whitelight".'

        defaults['method'] = "subpixel"
        dtypes['method'] = str
        descr['method'] = 'What method should be used to generate the datacube. There are currently two options: ' \
                          '(1) "subpixel" (default) - this algorithm divides each pixel in the spec2d frames ' \
                          'into subpixels, and assigns each subpixel to a voxel of the datacube. Flux is conserved, ' \
                          'but voxels are correlated, and the error spectrum does not account for covariance between ' \
                          'adjacent voxels. See also, spec_subpixel and spat_subpixel. ' \
                          '(2) "NGP" (nearest grid point) - this algorithm is effectively a 3D histogram. Flux is ' \
                          'conserved, voxels are not correlated, however this option suffers the same downsides as ' \
                          'any histogram; the choice of bin sizes can change how the datacube appears. This algorithm ' \
                          'takes each pixel on the spec2d frame and puts the flux of this pixel into one voxel in the ' \
                          'datacube. Depending on the binning used, some voxels may be empty (zero flux) while a ' \
                          'neighboring voxel might contain the flux from two spec2d pixels. Note that all spec2d ' \
                          'pixels that contribute to the same voxel are inverse variance weighted (e.g. if two ' \
                          'pixels have the same variance, the voxel would be assigned the average flux of the two ' \
                          'pixels).'
        # '(3) "resample" - this algorithm resamples the spec2d frames into a datacube. ' \
        # 'Flux is conserved, but voxels are correlated, and the error spectrum does not account ' \
        # 'for covariance between neighbouring pixels. ' \

        defaults['spec_subpixel'] = 5
        dtypes['spec_subpixel'] = int
        descr['spec_subpixel'] = 'When method=subpixel, spec_subpixel sets the subpixellation scale of ' \
                                 'each detector pixel in the spectral direction. The total number of subpixels ' \
                                 'in each pixel is given by spec_subpixel x spat_subpixel. The default option ' \
                                 'is to divide each spec2d pixel into 25 subpixels during datacube creation. ' \
                                 'See also, spat_subpixel.'

        defaults['spat_subpixel'] = 5
        dtypes['spat_subpixel'] = int
        descr['spat_subpixel'] = 'When method=subpixel, spat_subpixel sets the subpixellation scale of ' \
                                 'each detector pixel in the spatial direction. The total number of subpixels ' \
                                 'in each pixel is given by spec_subpixel x spat_subpixel. The default option ' \
                                 'is to divide each spec2d pixel into 25 subpixels during datacube creation. ' \
                                 'See also, spec_subpixel.'

        defaults['ra_min'] = None
        dtypes['ra_min'] = float
        descr['ra_min'] = 'Minimum RA to use when generating the WCS. If None, the default is minimum RA ' \
                          'based on the WCS of all spaxels. Units should be degrees.'

        defaults['ra_max'] = None
        dtypes['ra_max'] = float
        descr['ra_max'] = 'Maximum RA to use when generating the WCS. If None, the default is maximum RA ' \
                          'based on the WCS of all spaxels. Units should be degrees.'

        defaults['dec_min'] = None
        dtypes['dec_min'] = float
        descr['dec_min'] = 'Minimum DEC to use when generating the WCS. If None, the default is minimum DEC ' \
                           'based on the WCS of all spaxels. Units should be degrees.'

        defaults['dec_max'] = None
        dtypes['dec_max'] = float
        descr['dec_max'] = 'Maximum DEC to use when generating the WCS. If None, the default is maximum DEC ' \
                           'based on the WCS of all spaxels. Units should be degrees.'

        defaults['wave_min'] = None
        dtypes['wave_min'] = float
        descr['wave_min'] = 'Minimum wavelength to use when generating the WCS. If None, the default is ' \
                            'minimum wavelength based on the WCS of all spaxels. Units should be Angstroms.'

        defaults['wave_max'] = None
        dtypes['wave_max'] = float
        descr['wave_max'] = 'Maximum wavelength to use when generating the WCS. If None, the default is ' \
                            'maximum wavelength based on the WCS of all spaxels. Units should be Angstroms.'

        defaults['spatial_delta'] = None
        dtypes['spatial_delta'] = float
        descr['spatial_delta'] = 'The spatial size of each spaxel to use when generating the WCS (in arcsec). ' \
                                 'If None, the default is set by the spectrograph file.'

        defaults['wave_delta'] = None
        dtypes['wave_delta'] = float
        descr['wave_delta'] = 'The wavelength step to use when generating the WCS (in Angstroms). ' \
                                'If None, the default is set by the wavelength solution.'

        defaults['astrometric'] = True
        dtypes['astrometric'] = bool
        descr['astrometric'] = 'If true, an astrometric correction will be applied using the alignment frames.'

        defaults['grating_corr'] = True
        dtypes['grating_corr'] = bool
        descr['grating_corr'] = 'This option performs a small correction for the relative blaze function of all ' \
                                'input frames that have (even slightly) different grating angles, or if you are ' \
                                'flux calibrating your science data with a standard star that was observed with ' \
                                'a slightly different setup.'

        defaults['scale_corr'] = None
        dtypes['scale_corr'] = str
        descr['scale_corr'] = 'This option performs a small correction for the relative spectral illumination ' \
                              'scale of different spec2D files. Specify the relative path+file to the spec2D ' \
                              'file that you would like to use for the relative scaling. If you want to perform ' \
                              'this correction, it is best to use the spec2d file with the highest S/N sky spectrum. ' \
                              'You should choose the same frame for both the standards and science frames.'

        defaults['skysub_frame'] = 'image'
        dtypes['skysub_frame'] = str
        descr['skysub_frame'] = 'Set the sky subtraction to be implemented. The default behaviour is to subtract ' \
                                'the sky using the model that is derived from each individual image (i.e. set ' \
                                'this parameter to "image"). To turn off sky subtraction completely, set this ' \
                                'parameter to "none" (all lowercase). Finally, if you want to use a different frame ' \
                                'for the sky subtraction, specify the relative path+file to the spec2D file that you ' \
                                'would like to use for the sky subtraction. The model fit to the sky of the specified ' \
                                'frame will be used. Note, the sky and science frames do not need to have the same ' \
                                'exposure time.'

        # Instantiate the parameter set
        super(CubePar, self).__init__(list(pars.keys()),
                                      values=list(pars.values()),
                                      defaults=list(defaults.values()),
                                      options=list(options.values()),
                                      dtypes=list(dtypes.values()),
                                      descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])

        # Basic keywords
        parkeys = ['slit_spec', 'output_filename', 'standard_cube', 'reference_image', 'save_whitelight',
                   'method', 'spec_subpixel', 'spat_subpixel', 'ra_min', 'ra_max', 'dec_min', 'dec_max',
                   'wave_min', 'wave_max', 'spatial_delta', 'wave_delta', 'relative_weights', 'combine',
                   'astrometric', 'grating_corr', 'scale_corr', 'skysub_frame']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for CubePar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        allowed_methods = ["subpixel", "NGP"]#, "resample"
        if self.data['method'] not in allowed_methods:
            raise ValueError("The 'method' must be one of:\n"+", ".join(allowed_methods))


class FluxCalibratePar(ParSet):
    """
    A parameter set holding the arguments for how to perform the flux
    calibration.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    def __init__(self, extinct_correct=None, extinct_file=None, extrap_sens=None, use_archived_sens = False):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        defaults['extrap_sens'] = False
        dtypes['extrap_sens'] = bool
        descr['extrap_sens'] = "If False (default), the code will crash if one tries to use " \
                               "sensfunc at wavelengths outside its defined domain. By changing the " \
                               "par['sensfunc']['extrap_blu'] and par['sensfunc']['extrap_red'] this domain " \
                               "can be extended. If True the code will blindly extrapolate."


        defaults['extinct_correct'] = None
        dtypes['extinct_correct'] = bool
        descr['extinct_correct'] = 'The default behavior for atmospheric extinction corrections is that if UVIS algorithm is used ' \
                                   '(which does not correct for telluric absorption) than an atmospheric extinction model ' \
                                   'is used to correct for extinction below 10,000A, whereas if the IR algorithm is used, then ' \
                                   'no extinction correction is applied since the atmosphere is modeled directly. To follow these ' \
                                   'defaults based on the algorithm this parameter should be set to ``extinct_correct=None``. If instead this ' \
                                   'parameter is set, this overide this default behavior. In other words, it will force an extinction correction ' \
                                   'if ``extinct_correct=True``, and will not perform an extinction correction if ``extinct_correct=False``.'

        defaults['extinct_file'] = 'closest'
        dtypes['extinct_file'] = str
        descr['extinct_file'] = 'If ``extinct_file=\'closest\'`` the code will select the PypeIt-included extinction ' \
                                'file for the closest observatory (within 5 deg, geographic coordinates) to the telescope ' \
                                'identified in ``std_file`` (see :ref:`extinction_correction` for the list of currently '\
                                'included files).  If constructing a sesitivity function for a telescope not within 5 deg ' \
                                'of a listed observatory, this parameter may be set to the name of one of the listed ' \
                                'extinction files.  Alternatively, a custom extinction file may be installed in the ' \
                                'PypeIt cache using the ``pypeit_install_extinctfile`` script; this parameter may then ' \
                                'be set to the name of the custom extinction file.'

        defaults['use_archived_sens'] = False
        dtypes['use_archived_sens'] = bool
        descr['use_archived_sens'] = 'Use an archived sensfunc to flux calibration'

        # Instantiate the parameter set
        super(FluxCalibratePar, self).__init__(list(pars.keys()),
                                                 values=list(pars.values()),
                                                 defaults=list(defaults.values()),
                                                 dtypes=list(dtypes.values()),
                                                 descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])
        parkeys = ['extinct_correct', 'extinct_file', 'extrap_sens', 'use_archived_sens']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for FluxCalibratePar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        pass


class SensFuncPar(ParSet):
    """
    A parameter set holding the arguments for sensitivity function computation using the UV algorithm, see
    sensfunc.SensFuncUV

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    def __init__(self, extrap_blu=None, extrap_red=None, samp_fact=None, multi_spec_det=None, algorithm=None, UVIS=None,
                 IR=None, polyorder=None, star_type=None, star_mag=None, star_ra=None,
                 star_dec=None, mask_hydrogen_lines=None, mask_helium_lines=None, hydrogen_mask_wid=None):
        # Grab the parameter names and values from the function arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k, values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        defaults['extrap_blu'] = 0.1
        dtypes['extrap_blu'] = float
        descr['extrap_blu'] = 'Fraction of minimum wavelength coverage to grow the wavelength coverage of the ' \
                              'sensitivitity function in the blue direction (`i.e.`, if the standard star spectrum ' \
                              'cuts off at ``wave_min``) the sensfunc will be extrapolated to cover down to ' \
                              ' (1.0 - ``extrap_blu``) * ``wave_min``'


        defaults['extrap_red'] = 0.1
        dtypes['extrap_red'] = float
        descr['extrap_red'] = 'Fraction of maximum wavelength coverage to grow the wavelength coverage of the ' \
                              'sensitivitity function in the red direction (`i.e.`, if the standard star spectrum' \
                              'cuts off at ``wave_max``) the sensfunc will be extrapolated to cover up to ' \
                              ' (1.0 + ``extrap_red``) * ``wave_max``'

        defaults['samp_fact'] = 1.5
        dtypes['samp_fact'] = float
        descr['samp_fact'] = 'Sampling factor to make the wavelength grid for sensitivity function finer or coarser. ' \
                             'samp_fact > 1.0 oversamples (finer), samp_fact < 1.0 undersamples (coarser).'

        defaults['multi_spec_det'] = None
        dtypes['multi_spec_det'] = list
        descr['multi_spec_det'] = 'List of detectors (identified by their string name, like ' \
                                  'DET01) to splice together for multi-detector instruments ' \
                                  '(e.g. DEIMOS). It is assumed that there is *no* overlap in ' \
                                  'wavelength across detectors (might be ok if there is).  If ' \
                                  'entered as a list of integers, they should be converted to ' \
                                  'the detector name.  **Cannot be used with detector mosaics.**'

        defaults['algorithm'] = 'UVIS'
        dtypes['algorithm'] = str
        options['algorithm'] = SensFuncPar.valid_algorithms()
        descr['algorithm'] = "Specify the algorithm for computing the sensitivity function. The options are: " \
                             r" (1) UVIS = Should be used for data with :math:`\lambda < 7000` A. " \
                             "No detailed model of telluric absorption but corrects for atmospheric extinction. " \
                             r" (2) IR = Should be used for data with :math:`\lambda > 7000` A. " \
                             "Peforms joint fit for sensitivity function and telluric absorption using HITRAN models."


        defaults['UVIS'] = SensfuncUVISPar()
        dtypes['UVIS'] = [ParSet, dict ]
        descr['UVIS'] = 'Parameters for the UVIS sensfunc algorithm'

        defaults['IR'] = TelluricPar()
        dtypes['IR'] = [ ParSet, dict ]
        descr['IR'] = 'Parameters for the IR sensfunc algorithm'

        # JFH SHould the default by higher like 8?
        defaults['polyorder'] = 5
        dtypes['polyorder'] = [int, list]
        descr['polyorder'] = 'Polynomial order for sensitivity function fitting'

        defaults['star_type'] = None
        dtypes['star_type'] = str
        descr['star_type'] = 'Spectral type of the standard star (for near-IR mainly)'

        defaults['star_mag'] = None
        dtypes['star_mag'] = float
        descr['star_mag'] = 'Magnitude of the standard star (for near-IR mainly)'

        defaults['star_ra'] = None
        dtypes['star_ra'] = float
        descr['star_ra'] = 'RA of the standard star. This will override values in the header (`i.e.`, if they are wrong or absent)'

        defaults['star_dec'] = None
        dtypes['star_dec'] = float
        descr['star_dec'] = 'DEC of the standard star. This will override values in the header (`i.e.`, if they are wrong or absent)'

        defaults['mask_hydrogen_lines'] = True
        dtypes['mask_hydrogen_lines'] = bool
        descr['mask_hydrogen_lines'] = 'Mask hydrogen Balmer, Paschen, Brackett, and Pfund recombination lines in the sensitivity function fit. ' \
                                       'A region equal to ``hydrogen_mask_wid`` on either side of the line center is masked.'

        defaults['hydrogen_mask_wid'] = 10.0
        dtypes['hydrogen_mask_wid'] = float
        descr['hydrogen_mask_wid'] = 'Mask width from line center for hydrogen recombination lines in Angstroms (total mask width is 2x this value).'

        defaults['mask_helium_lines'] = False
        dtypes['mask_helium_lines'] = bool
        descr['mask_helium_lines'] = 'Mask certain ``HeII`` recombination lines prominent in O-type stars in the sensitivity function fit ' \
                                     'A region equal to 0.5 * ``hydrogen_mask_wid`` on either side of the line center is masked.'

        # Instantiate the parameter set
        super(SensFuncPar, self).__init__(list(pars.keys()),
                                          values=list(pars.values()),
                                          defaults=list(defaults.values()),
                                          options=list(options.values()),
                                          dtypes=list(dtypes.values()),
                                          descr=list(descr.values()))
#        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])

        # Single element parameters
        parkeys = ['extrap_blu', 'extrap_red', 'samp_fact', 'multi_spec_det', 'algorithm',
                   'polyorder', 'star_type', 'star_mag', 'star_ra', 'star_dec',
                   'mask_hydrogen_lines', 'mask_helium_lines', 'hydrogen_mask_wid']

        # All parameters, including nested ParSets
        allkeys = parkeys + ['UVIS', 'IR']

        badkeys = np.array([pk not in allkeys for pk in k])
        if np.any(badkeys):
            raise ValueError(f"{','.join(k[badkeys])} are not recognized key(s) for SensFuncPar.")

        kwargs = {}
        # Single element parameters
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        # Parameters that are themselves ParSets
        pk = 'UVIS'
        kwargs[pk] = SensfuncUVISPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'IR'
        kwargs[pk] = TelluricPar.from_dict(cfg[pk]) if pk in k else None

        return cls(**kwargs)

    @staticmethod
    def valid_algorithms():
        """
        Return the valid sensitivity algorithms.
        """
        return ['UVIS', 'IR']


class SensfuncUVISPar(ParSet):
    """
    A parameter set holding the arguments for sensitivity function computation using the UV algorithm, see
    sensfunc.SensFuncUV

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    def __init__(self, std_file=None, std_obj_id=None, sensfunc=None, extinct_correct=None,
                 extinct_file=None, telluric_correct=None, telluric=None, polycorrect=None,
                 polyfunc=None, nresln=None, resolution=None, trans_thresh=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())


        # These are the UV sensfunc parameters
        dtypes['std_file'] = str
        descr['std_file'] = 'Standard star file to generate sensfunc'

        dtypes['std_obj_id'] = [str, int]
        descr['std_obj_id'] = 'Specifies object in spec1d file to use as standard.' \
            ' The brightest object found is used otherwise.'


        dtypes['sensfunc'] = str
        descr['sensfunc'] = 'FITS file that contains or will contain the sensitivity function.'


        defaults['extinct_correct'] = True
        dtypes['extinct_correct'] = bool
        descr['extinct_correct'] = 'If ``extinct_correct=True`` the code will use an atmospheric extinction model to ' \
                                   'extinction correct the data below 10000A. Note that this correction makes no ' \
                                   'sense if one is telluric correcting and this shold be set to False'

        defaults['extinct_file'] = 'closest'
        dtypes['extinct_file'] = str
        descr['extinct_file'] = 'If ``extinct_file=\'closest\'`` the code will select the PypeIt-included extinction ' \
                                'file for the closest observatory (within 5 deg, geographic coordinates) to the telescope ' \
                                'identified in ``std_file`` (see :ref:`extinction_correction` for the list of currently '\
                                'included files).  If constructing a sesitivity function for a telescope not within 5 deg ' \
                                'of a listed observatory, this parameter may be set to the name of one of the listed ' \
                                'extinction files.  Alternatively, a custom extinction file may be installed in the ' \
                                'PypeIt cache using the ``pypeit_install_extinctfile`` script; this parameter may then ' \
                                'be set to the name of the custom extinction file.'

        defaults['telluric_correct'] = False
        dtypes['telluric_correct'] = bool
        descr['telluric_correct'] = "If ``telluric_correct=True`` the code will grab the sens_dict['telluric'] tag from the " \
                                    "sensfunc dictionary and apply it to the data."

        defaults['telluric'] = False
        dtypes['telluric'] = bool
        descr['telluric'] = 'If ``telluric=True`` the code creates a synthetic standard star spectrum using the Kurucz models, ' \
            'the sens func is created setting nresln=1.5 it contains the correction for telluric lines.'

        defaults['polycorrect'] = True
        dtypes['polycorrect'] = bool
        descr['polycorrect'] = 'Whether you want to correct the sensfunc with polynomial in the telluric and recombination line regions'

        defaults['polyfunc'] = False
        dtypes['polyfunc'] = bool
        descr['polyfunc'] = 'Whether you want to use the polynomial fit as your final SENSFUNC'


        defaults['nresln'] = 20
        dtypes['nresln'] = [int, float]
        descr['nresln'] = 'Parameter governing the spacing of the bspline breakpoints in terms of number of resolution elements.'


        defaults['resolution'] = 3000.0
        dtypes['resolution'] = [int, float]
        descr['resolution'] = 'Expected resolution of the standard star spectrum. This should be measured from the data.'

        defaults['trans_thresh'] = 0.9
        dtypes['trans_thresh'] = float
        descr['trans_thresh'] = 'Parameter for selecting telluric regions which are masked. Locations below this ' \
                                'transmission value are masked. If you have significant telluric absorption you should ' \
                                'be using telluric.sensnfunc_telluric'

        # Instantiate the parameter set
        super(SensfuncUVISPar, self).__init__(list(pars.keys()),
                                          values=list(pars.values()),
                                          defaults=list(defaults.values()),
                                          dtypes=list(dtypes.values()),
                                          descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])
        parkeys = ['sensfunc', 'extinct_correct', 'extinct_file', 'telluric_correct', 'std_file',
                   'std_obj_id', 'telluric', 'polyfunc','polycorrect', 'nresln', 'resolution', 'trans_thresh']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for SensfuncUVISPar.'.format(k[badkeys]))

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

class SlitMaskPar(ParSet):
    """
    A parameter set holding the arguments for fussing with
    slitmask ingestion and object assignment

    A list of these objects can be included in an instance of
    :class:`SlitMaskPar` to perform a set of user-defined
    extractions.

    Args:


    """
    def __init__(self, obj_toler=None, assign_obj=None, snr_thrshd=None,
                 slitmask_offset=None, use_dither_offset=None, bright_maskdef_id=None, extract_missing_objs=None,
                 missing_objs_fwhm=None, missing_objs_boxcar_rad=None, use_alignbox=None):

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

        defaults['obj_toler'] = 1.
        dtypes['obj_toler'] = [int, float]
        descr['obj_toler'] = 'If slitmask design information is provided, and slit matching is performed ' \
                             '(``use_maskdesign = True`` in ``EdgeTracePar``), this parameter provides ' \
                             'the desired tolerance (arcsec) to match sources to targeted objects'

        defaults['assign_obj'] = False
        dtypes['assign_obj'] = bool
        descr['assign_obj'] = 'If SlitMask object was generated, assign RA,DEC,name to detected objects'

        defaults['use_alignbox'] = False
        dtypes['use_alignbox'] = bool
        descr['use_alignbox'] = 'Use stars in alignment boxes to compute the slitmask offset. ' \
                                'If this is set to ``True`` PypeIt will NOT compute ' \
                                'the offset using ``snr_thrshd`` or ``bright_maskdef_id``'

        defaults['snr_thrshd'] = 50.
        dtypes['snr_thrshd'] = [int, float]
        descr['snr_thrshd'] = 'Objects detected above this S/N threshold will ' \
                               'be used to compute the slitmask offset. This is the default behaviour for DEIMOS ' \
                               ' unless ``slitmask_offset``, ``bright_maskdef_id`` or ``use_alignbox`` is set.'

        defaults['slitmask_offset'] = None
        dtypes['slitmask_offset'] = [int, float]
        descr['slitmask_offset'] = 'User-provided slitmask offset (pixels) from the position expected by ' \
                                   'the slitmask design. This is optional, and if set PypeIt will NOT compute ' \
                                   'the offset using ``snr_thrshd`` or ``bright_maskdef_id``.'

        defaults['use_dither_offset'] = False
        dtypes['use_dither_offset'] = bool
        descr['use_dither_offset'] = 'Use the dither offset recorded in the header of science frames as the value ' \
                                     'of the slitmask offset. This is currently only available for Keck MOSFIRE ' \
                                     'reduction and it is set as the default for this instrument. If set PypeIt will ' \
                                     'NOT compute the offset using ``snr_thrshd`` or ``bright_maskdef_id``. ' \
                                     'However, it is ignored if ``slitmask_offset`` is provided. '

        defaults['bright_maskdef_id'] = None
        dtypes['bright_maskdef_id'] = int
        descr['bright_maskdef_id'] = '`maskdef_id` (corresponding to `dSlitId` and `Slit_Number` in the DEIMOS ' \
                                     'and MOSFIRE slitmask design, respectively) of a ' \
                                     'slit containing a bright object that will be used to compute the ' \
                                     'slitmask offset. This parameter is optional and is ignored ' \
                                     'if ``slitmask_offset`` is provided.'

        defaults['extract_missing_objs'] = False
        dtypes['extract_missing_objs'] = bool
        descr['extract_missing_objs'] = 'Force extraction of undetected objects at the location expected ' \
                                        'from the slitmask design.'

        defaults['missing_objs_fwhm'] = None
        dtypes['missing_objs_fwhm'] = [int, float]
        descr['missing_objs_fwhm'] = 'Indicates the FWHM in arcsec for the force extraction of undetected objects. ' \
                                     'PypeIt will try to determine the FWHM from the flux profile ' \
                                     '(by using ``missing_objs_fwhm`` as initial guess). ' \
                                     'If the FWHM cannot be determined, ``missing_objs_fwhm`` will be assumed. ' \
                                     'If you do not want PypeIt to try to determine the FWHM set the ' \
                                     'parameter ``use_user_fwhm`` in ``ExtractionPar`` to True. ' \
                                     'If ``missing_objs_fwhm`` is ``None`` (which is the default) PypeIt will use ' \
                                     'the median FWHM of all the detected objects.'

        defaults['missing_objs_boxcar_rad'] = 1.0
        dtypes['missing_objs_boxcar_rad'] = [int, float]
        descr['missing_objs_boxcar_rad'] = 'Indicates the boxcar radius in arcsec for the force ' \
                                           'extraction of undetected objects. '

        # Instantiate the parameter set
        super(SlitMaskPar, self).__init__(list(pars.keys()),
                                          values=list(pars.values()),
                                          defaults=list(defaults.values()),
                                          dtypes=list(dtypes.values()),
                                          descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])
        parkeys = ['obj_toler', 'assign_obj', 'snr_thrshd', 'slitmask_offset', 'use_dither_offset',
                   'bright_maskdef_id', 'extract_missing_objs', 'missing_objs_fwhm', 'missing_objs_boxcar_rad',
                   'use_alignbox']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for SlitMaskPar.'.format(
                                k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        pass


class TelluricPar(ParSet):
    """
    A parameter set holding the arguments for sensitivity function computation using the UV algorithm, see
    sensfunc.SensFuncUV

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """

    def __init__(self, telgridfile=None, sn_clip=None, resln_guess=None, resln_frac_bounds=None, pix_shift_bounds=None,
                 delta_coeff_bounds=None, minmax_coeff_bounds=None, maxiter=None,
                 sticky=None, lower=None, upper=None, seed=None, tol=None, popsize=None, recombination=None, polish=None,
                 disp=None, objmodel=None, redshift=None, delta_redshift=None, pca_file=None, npca=None,
                 bal_wv_min_max=None, bounds_norm=None, tell_norm_thresh=None, only_orders=None, pca_lower=None,
                 pca_upper=None, star_type=None, star_mag=None, star_ra=None, star_dec=None,
                 func=None, model=None, polyorder=None, fit_wv_min_max=None, mask_lyman_a=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k, values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        defaults['telgridfile'] = None
        dtypes['telgridfile'] = str
        descr['telgridfile'] = 'File containing the telluric grid for the observatory in question. These grids are ' \
                               'generated from HITRAN models for each observatory using nominal site parameters. They ' \
                               'must be downloaded from the GoogleDrive and installed in your PypeIt installation via ' \
                               'the pypeit_install_telluric script. NOTE: This parameter no longer includes the full ' \
                               'pathname to the Telluric Grid file, but is just the filename of the grid itself.'

        defaults['sn_clip'] = 30.0
        dtypes['sn_clip'] = [int, float]
        descr['sn_clip'] = 'This adds an error floor to the ivar, preventing too much rejection at high-S/N (`i.e.`, ' \
                          'standard stars, bright objects) using the function utils.clip_ivar. A small erorr is added ' \
                          'to the input ivar so that the output ivar_out will never give S/N greater than sn_clip. This ' \
                          'prevents overly aggressive rejection in high S/N ratio spectra which neverthless differ at a ' \
                          'level greater than the formal S/N due to the fact that our telluric models are only good to ' \
                          'about 3%.'

        defaults['resln_guess'] = None
        dtypes['resln_guess'] = [int, float]
        descr['resln_guess'] = 'A guess for the resolution of your spectrum expressed as lambda/dlambda. The resolution ' \
                               'is fit explicitly as part of the telluric model fitting, but this guess helps determine ' \
                               'the bounds for the optimization (see next). If not provided, the  wavelength sampling of ' \
                               'your spectrum will be used and the resolution calculated using a typical sampling of 3 ' \
                               'spectral pixels per resolution element.'


        pars['resln_frac_bounds'] = tuple_force(pars['resln_frac_bounds'])
        defaults['resln_frac_bounds'] = (0.5,1.5)
        dtypes['resln_frac_bounds'] = tuple
        descr['resln_frac_bounds'] = 'Bounds for the resolution fit optimization which is part of the telluric model. ' \
                                     'This range is in units of the resln_guess, so the (0.5, 1.5) would bound the ' \
                                     'spectral resolution fit to be within the range ' \
                                     'bounds_resln = (0.5*resln_guess, 1.5*resln_guess)'

        pars['pix_shift_bounds'] = tuple_force(pars['pix_shift_bounds'])
        defaults['pix_shift_bounds'] = (-5.0,5.0)
        dtypes['pix_shift_bounds'] = tuple
        descr['pix_shift_bounds'] = 'Bounds for the pixel shift optimization in telluric model fit in units of pixels. ' \
                                    'The atmosphere will be allowed to shift within this range during the fit.'

        pars['delta_coeff_bounds'] = tuple_force(pars['delta_coeff_bounds'])
        defaults['delta_coeff_bounds'] = (-20.0, 20.0)
        dtypes['delta_coeff_bounds'] = tuple
        descr['delta_coeff_bounds'] = 'Parameters setting the polynomial coefficient bounds for sensfunc optimization.'

        pars['minmax_coeff_bounds'] = tuple_force(pars['minmax_coeff_bounds'])
        defaults['minmax_coeff_bounds'] = (-5.0, 5.0)
        dtypes['minmax_coeff_bounds'] = tuple
        descr['minmax_coeff_bounds'] = "Parameters setting the polynomial coefficient bounds for sensfunc optimization. " \
                                       "Bounds are currently determined as follows. We compute an initial fit to the " \
                                       "sensfunc in the pypeit.core.telluric.init_sensfunc_model function. That deterines " \
                                       "a set of coefficients. The bounds are then determined according to: " \
                                       "[(np.fmin(np.abs(this_coeff)*obj_params['delta_coeff_bounds'][0], " \
                                       "obj_params['minmax_coeff_bounds'][0]), " \
                                       "np.fmax(np.abs(this_coeff)*obj_params['delta_coeff_bounds'][1], " \
                                       "obj_params['minmax_coeff_bounds'][1]))]"

        defaults['maxiter'] = 2
        dtypes['maxiter'] = int
        descr['maxiter'] = 'Maximum number of iterations for the telluric + object model fitting. The code performs ' \
                           'multiple iterations rejecting outliers at each step. The fit is then performed anew to the ' \
                           'remaining good pixels. For this reason if you run with the disp=True option, you will see ' \
                           'that the f(x) loss function gets progressively better during the iterations.'

        defaults['sticky'] = True
        dtypes['sticky'] = bool
        descr['sticky'] = 'Sticky parameter for the utils.djs_reject algorithm for iterative model fit rejection.  ' \
                          'If set to True then points rejected from a previous iteration are kept rejected, in other ' \
                          'words the bad pixel mask is the OR of all previous iterations and rejected pixels accumulate. ' \
                          'If set to False, the bad pixel mask is the mask from the previous iteration, and if the model ' \
                          'fit changes between iterations, points can alternate from being rejected to not rejected. ' \
                          'At present this code only performs optimizations with differential evolution and experience ' \
                          'shows that sticky needs to be True in order for these to converge. This is because the ' \
                          'outliers can be so large that they dominate the loss function, and one never iteratively ' \
                          'converges to a good model fit. In other words, the deformations in the model between ' \
                          'iterations with sticky=False are too small to approach a reasonable fit.'

        defaults['lower'] = 3.0
        dtypes['lower'] = [int, float]
        descr['lower'] = 'Lower rejection threshold in units of sigma_corr*sigma, where sigma is the formal noise of the ' \
                         'spectrum, and sigma_corr is an empirically determined correction to the formal error. The ' \
                         'distribution of input chi (defined by chi = (data - model)/sigma) values is analyzed, and a ' \
                         'correction factor to the formal error sigma_corr is returned which is multiplied into the ' \
                         'formal errors. In this way, a rejection threshold of i.e. 3-sigma, will always correspond to ' \
                         'roughly the same percentile.  This renormalization is performed with ' \
                         'coadd1d.renormalize_errors function, and guarantees that rejection is not too agressive in ' \
                         'cases where the empirical errors determined from the chi-distribution differ significantly ' \
                         'from the formal noise which is used to determine chi.'

        defaults['upper'] = 3.0
        dtypes['upper'] = [int, float]
        descr['upper'] = 'Upper rejection threshold in units of sigma_corr*sigma, where sigma is the formal noise of the ' \
                         'spectrum, and sigma_corr is an empirically determined correction to the formal error. See ' \
                         'above for description.'

        defaults['seed'] = 777
        dtypes['seed'] = int
        descr['seed'] = 'An initial seed for the differential evolution optimization, which is a random process. ' \
                        'The default is a seed = 777 which will be used to generate a unique seed for every order. ' \
                        'A specific seed is used because otherwise the random number generator will use the time for ' \
                        'the seed, and the results will not be reproducible.'


        defaults['tol'] = 1e-3
        dtypes['tol'] = float
        descr['tol'] = 'Relative tolerance for converage of the differential evolution optimization. See ' \
                       'scipy.optimize.differential_evolution for details.'


        defaults['popsize'] = 30
        dtypes['popsize'] = int
        descr['popsize'] = 'A multiplier for setting the total population size for the differential evolution ' \
                           'optimization. See scipy.optimize.differential_evolution for details.'

        defaults['recombination'] = 0.7
        dtypes['recombination'] = [int, float]
        descr['recombination'] = 'The recombination constant for the differential evolution optimization. This should ' \
                                 'be in the range [0, 1]. See scipy.optimize.differential_evolution for details.'

        defaults['polish'] = True
        dtypes['polish'] = bool
        descr['polish'] = 'If True then differential evolution will perform an additional optimizatino at the end to ' \
                          'polish the best fit at the end, which can improve the optimization slightly. See ' \
                          'scipy.optimize.differential_evolution for details.'

        defaults['disp'] = False
        dtypes['disp'] = bool
        descr['disp'] = 'Argument for scipy.optimize.differential_evolution which will display status messages to the ' \
                        'screen indicating the status of the optimization. See documentation for telluric.Telluric ' \
                        'for a description of the output and how to know if things are working well.'


        defaults['only_orders'] = None
        dtypes['only_orders'] = [int, list, np.ndarray]
        descr['only_orders'] = "Order number, or list of order numbers if you only want to fit specific orders"

        defaults['objmodel'] = None
        dtypes['objmodel'] = str
        descr['objmodel'] = 'The object model to be used for telluric fitting. Currently the options are: qso, star, and poly'

        ### Parameters for qso_telluric
        defaults['redshift'] = 0.0
        dtypes['redshift'] = [int, float]
        descr['redshift'] = 'The redshift for the object model. This is currently only used by objmodel=qso'

        defaults['delta_redshift'] = 0.1
        dtypes['delta_redshift'] = float
        descr['delta_redshift'] = 'Range within the redshift can be varied for telluric fitting, i.e. the code performs a bounded optimization within ' \
                                  'the redshift +- delta_redshift'

        defaults['pca_file'] = 'qso_pca_1200_3100.fits'
        dtypes['pca_file'] = str
        descr['pca_file'] = 'Fits file containing quasar PCA model. Needed for objmodel=qso.  NOTE: This parameter no longer includes the full ' \
                               'pathname to the Telluric Model file, but is just the filename of the model itself.'

        defaults['npca'] = 8
        dtypes['npca'] = int
        descr['npca'] = 'Number of pca for the objmodel=qso qso PCA fit'

        defaults['bal_wv_min_max'] = None
        dtypes['bal_wv_min_max'] = [list, np.ndarray]
        descr['bal_wv_min_max'] = 'Min/max wavelength of broad absorption features. If there are several BAL features, ' \
                            'the format for this mask is [wave_min_bal1, wave_max_bal1,wave_min_bal2, ' \
                            'wave_max_bal2,...]. These masked pixels will be ignored during the fitting.'

        pars['bounds_norm'] = tuple_force(pars['bounds_norm'])
        defaults['bounds_norm'] = (0.1, 3.0)
        dtypes['bounds_norm'] = tuple
        descr['bounds_norm'] = "Normalization bounds for scaling the initial object model."

        defaults['tell_norm_thresh'] = 0.9
        dtypes['tell_norm_thresh'] = [int, float]
        descr['tell_norm_thresh'] = "Threshold of telluric absorption region"

        defaults['pca_lower'] = 1220.0
        dtypes['pca_lower'] = [int, float]
        descr['pca_lower'] = "Minimum wavelength for the qso pca model"

        defaults['pca_upper'] = 3100.0
        dtypes['pca_upper'] = [int, float]
        descr['pca_upper'] = "Maximum wavelength for the qso pca model"

        defaults['mask_lyman_a'] = True
        dtypes['mask_lyman_a'] = bool
        descr['mask_lyman_a'] = 'Mask the blueward of Lyman-alpha line during the fitting?'

        ### Start parameters for star_telluric
        defaults['star_type'] = None
        dtypes['star_type'] = str
        descr['star_type'] = 'stellar type'

        defaults['star_mag'] = None
        dtypes['star_mag'] = [float, int]
        descr['star_mag'] = 'AB magnitude in V band'

        defaults['star_ra'] = None
        dtypes['star_ra'] = float
        descr['star_ra'] = 'Object right-ascension in decimal deg'

        defaults['star_dec'] = None
        dtypes['star_dec'] = float
        descr['star_dec'] = 'Object declination in decimal deg'

        ### parameters for both star_telluric and poly_telluric
        defaults['func'] = 'legendre'
        dtypes['func'] = str
        descr['func'] = 'Polynomial model function'

        defaults['model'] = 'exp'
        dtypes['model'] = str
        descr['model'] = 'Types of polynomial model. Options are poly, square, exp corresponding to normal polynomial, '\
                         'squared polynomial, or exponentiated polynomial'

        defaults['polyorder'] = 3
        dtypes['polyorder'] = int
        descr['polyorder'] = "Order of the polynomial model fit"

        ### Start parameters for poly_telluric
        defaults['fit_wv_min_max'] = None
        dtypes['fit_wv_min_max'] = list
        descr['fit_wv_min_max'] = "Pixels within this mask will be used during the fitting. The format "\
                                   "is the same with bal_wv_min_max, but this mask is good pixel masks."


        # Instantiate the parameter set
        super(TelluricPar, self).__init__(list(pars.keys()),
                                          values=list(pars.values()),
                                          defaults=list(defaults.values()),
                                          dtypes=list(dtypes.values()),
                                          descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])
        parkeys = ['telgridfile', 'sn_clip', 'resln_guess', 'resln_frac_bounds',
                   'pix_shift_bounds', 'delta_coeff_bounds', 'minmax_coeff_bounds',
                   'maxiter', 'sticky', 'lower', 'upper', 'seed', 'tol',
                   'popsize', 'recombination', 'polish', 'disp', 'objmodel','redshift', 'delta_redshift',
                   'pca_file', 'npca', 'bal_wv_min_max', 'bounds_norm',
                   'tell_norm_thresh', 'only_orders', 'pca_lower', 'pca_upper',
                   'star_type','star_mag','star_ra','star_dec',
                   'func','model','polyorder','fit_wv_min_max','mask_lyman_a']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for TelluricPar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        pass
        # JFH add something in here which checks that the recombination value provided is bewteen 0 and 1, although
        # scipy.optimize.differential_evoluiton probalby checks this.


class ReduxPar(ParSet):
    """
    The parameter set used to hold arguments for functionality relevant
    to the overal reduction of the the data.

    Critically, this parameter set defines the spectrograph that was
    used to collect the data and the overall pipeline used in the
    reductions.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    def __init__(self, spectrograph=None, detnum=None, sortroot=None, calwin=None, scidir=None,
                 qadir=None, redux_path=None, ignore_bad_headers=None, slitspatnum=None,
                 maskIDs=None, quicklook=None):

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

        # NOTE: The validity of the spectrograph is checked by
        # load_spectrograph, so the specification of the viable options here is
        # not really necessary.
#        options['spectrograph'] = ReduxPar.valid_spectrographs()
        dtypes['spectrograph'] = str
        descr['spectrograph'] = 'Spectrograph that provided the data to be reduced.  ' \
                                'See :ref:`instruments` for valid options.'
#                                'Options are: {0}'.format(', '.join(options['spectrograph']))

        defaults['quicklook'] = False
        dtypes['quicklook'] = bool
        descr['quicklook'] = 'Run a quick look reduction? This is usually good if you want to quickly ' \
                             'reduce the data (usually at the telescope in real time) to get an initial ' \
                             'estimate of the data quality.'

        dtypes['detnum'] = [int, list]
        descr['detnum'] = 'Restrict reduction to a list of detector indices. ' \
                          'In case of mosaic reduction (currently only available for ' \
                          'Gemini/GMOS and Keck/DEIMOS) ``detnum`` should be a list of ' \
                          'tuples of the detector indices that are mosaiced together. ' \
                          'E.g., for Gemini/GMOS ``detnum`` would be ``[(1,2,3)]`` and for ' \
                          'Keck/DEIMOS it would be ``[(1, 5), (2, 6), (3, 7), (4, 8)]``' 

        dtypes['slitspatnum'] = [str, list]
        descr['slitspatnum'] = 'Restrict reduction to a set of slit DET:SPAT values (closest slit is used). ' \
                               'Example syntax -- slitspatnum = DET01:175,DET01:205 or MSC02:2234  If you are re-running the code, ' \
                               '(i.e. modifying one slit) you *must* have the precise SPAT_ID index.' 

        dtypes['maskIDs'] = [str, int, list]
        descr['maskIDs'] = 'Restrict reduction to a set of slitmask IDs ' \
                               'Example syntax -- ``maskIDs = 818006,818015`` ' \
                               'This must be used with detnum (for now).'

        dtypes['sortroot'] = str
        descr['sortroot'] = 'A filename given to output the details of the sorted files.  If ' \
                            'None, the default is the root name of the pypeit file.  If off, ' \
                            'no output is produced.'

        # TODO: Allow this to apply to each calibration frame type
        defaults['calwin'] = 0
        dtypes['calwin']   = [int, float]
        descr['calwin'] = 'The window of time in hours to search for calibration frames for a ' \
                          'science frame'

        # TODO: Explain what this actually does in the description.
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
        descr['redux_path'] = 'Path to folder for performing reductions.  Default is the ' \
                              'current working directory.'

        # Instantiate the parameter set
        super(ReduxPar, self).__init__(list(pars.keys()),
                                        values=list(pars.values()),
                                        defaults=list(defaults.values()),
                                        options=list(options.values()),
                                        dtypes=list(dtypes.values()),
                                        descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])

        # Basic keywords
        parkeys = [ 'spectrograph', 'quicklook', 'detnum', 'sortroot', 'calwin', 'scidir', 'qadir',
                    'redux_path', 'ignore_bad_headers', 'slitspatnum', 'maskIDs']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for ReduxPar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        # Finish
        return cls(**kwargs)

    def validate(self):
        if self.data['slitspatnum'] is not None:
            if self.data['maskIDs'] is not None:
                raise ValueError("You cannot assign both splitspatnum and maskIDs")
        if self.data['maskIDs'] is not None:
            if self.data['detnum'] is None:
                raise ValueError("You must assign detnum with maskIDs (for now)")
            # Recast as a list
            if not isinstance(self.data['maskIDs'], list):
                self.data['maskIDs'] = [self.data['maskIDs']]



class WavelengthSolutionPar(ParSet):
    """
    The parameter set used to hold arguments for the determination of
    wavelength solution.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    def __init__(self, reference=None, method=None, echelle=None, ech_nspec_coeff=None, ech_norder_coeff=None, ech_sigrej=None, lamps=None,
                 sigdetect=None, fwhm=None, fwhm_fromlines=None, reid_arxiv=None,
                 nreid_min=None, cc_thresh=None, cc_local_thresh=None, nlocal_cc=None,
                 rms_threshold=None, match_toler=None, func=None, n_first=None, n_final=None,
                 sigrej_first=None, sigrej_final=None, numsearch=None,
                 nfitpix=None, refframe=None,
                 nsnippet=None, use_instr_flag=None, wvrng_arxiv=None,
                 ech_separate_2d=None):

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
        descr['method'] = 'Method to use to fit the individual arc lines.  Note that some of ' \
                          'the available methods should not be used; they are unstable and ' \
                          'require significant parameter tweaking to succeed.  You should use ' \
                          'one of \'holy-grail\', \'reidentify\', or \'full_template\'.  ' \
                          '\'holy-grail\' attempts to get a first guess at line IDs by looking ' \
                          'for patterns in the line locations.  It is fully automated.  When ' \
                          'it works, it works well; however, it can fail catastrophically.  ' \
                          'Instead, \'reidentify\' and \'full_template\' are the preferred ' \
                          'methods.  They require an archived wavelength solution for your ' \
                          'specific instrument/grating combination as a reference.  ' \
                          'This is used to anchor the wavelength solution for the data being ' \
                          f"reduced.  All options are: {', '.join(options['method'])}."

        # Echelle wavelength calibration stuff
        # TODO: Is this needed? I.e., where do we need this parameter
        # when we don't have access to spectrograph.pypeline?
        defaults['echelle'] = False
        dtypes['echelle'] = bool
        descr['echelle'] = 'Is this an echelle spectrograph? If yes an additional 2-d fit ' \
                           'wavelength fit will be performed as a function of spectral pixel ' \
                           'and order number to improve the wavelength solution'

        defaults['ech_nspec_coeff'] = 4
        dtypes['ech_nspec_coeff'] = int
        descr['ech_nspec_coeff'] = 'For echelle spectrographs, this is the order of the final ' \
                                   '2d fit to the spectral dimension.  You should choose this ' \
                                   'to be the n_final of the fits to the individual orders.'

        defaults['ech_norder_coeff'] = 4
        dtypes['ech_norder_coeff'] = int
        descr['ech_norder_coeff'] = 'For echelle spectrographs, this is the order of the final ' \
                                    '2d fit to the order dimension.'

        defaults['ech_sigrej'] = 2.0
        dtypes['ech_sigrej'] = [int,float]
        descr['ech_sigrej'] = 'For echelle spectrographs, this is the sigma-clipping rejection ' \
                              'threshold in the 2d fit to spectral and order dimensions'

        defaults['ech_separate_2d'] = False
        dtypes['ech_separate_2d'] = bool
        descr['ech_separate_2d'] = 'For echelle spectrographs, fit the 2D solutions on separate detectors separately'

        # TODO: These needs to be tidied up so we can check for valid
        # lamps. Right now I'm not checking.

        # Force lamps to be a list
        if pars['lamps'] is not None and not isinstance(pars['lamps'], list):
            pars['lamps'] = [pars['lamps']]
        options['lamps'] = None
        #options['lamps'] = WavelengthSolutionPar.valid_lamps()
        dtypes['lamps'] = list
        descr['lamps'] = 'Name of one or more ions used for the wavelength calibration.  Use ' \
                         '``None`` for no calibration. Choose ``use_header`` to use the list of lamps ' \
                         'recorded in the header of the arc frames (this is currently ' \
                         'available only for Keck DEIMOS and LDT DeVeny).' # \
#                         'Options are: {0}'.format(', '.join(WavelengthSolutionPar.valid_lamps()))

        defaults['use_instr_flag'] = False
        dtypes['use_instr_flag'] = bool
        descr['use_instr_flag'] = 'If True, restrict to lines matching the instrument.  WARNING: This ' \
            'is only implemented for shane_kast_red + HolyGrail.  Do not use it unless you really know what you are doing.'

        # ToDo Should this be in counts or ADU? Currently the arcs are in ADU (which actually sort of makes sense here) but the
        # name of the parameter is counts. Perhaps we should just change this to nonlinear_adu or something to avoid confusion.

        # These are the parameters used for arc line detection
        # TODO: Why is this not always defined by the detectors of the
        # spectrograph?
        #defaults['nonlinear_counts'] = None
        #dtypes['nonlinear_counts'] = float
        #descr['nonlinear_counts'] = 'Arc lines above this saturation threshold are not used in wavelength solution fits because they cannot' \
        #                            'be accurately centroided'

        defaults['sigdetect'] = 5.
        dtypes['sigdetect'] =  [int, float, list, np.ndarray]
        descr['sigdetect'] = 'Sigma threshold above fluctuations for arc-line detection.  Arcs ' \
                             'are continuum subtracted and the fluctuations are computed after ' \
                             'continuum subtraction.  This can be a single number or a vector ' \
                             '(list or numpy array) that provides the detection threshold for ' \
                             'each slit.'

        defaults['fwhm'] = 4.
        dtypes['fwhm'] = [int, float]
        descr['fwhm'] = 'Spectral sampling of the arc lines. This is the FWHM of an arcline in ' \
                        'binned pixels of the input arc image'

        defaults['fwhm_fromlines'] = False
        dtypes['fwhm_fromlines'] = bool
        descr['fwhm_fromlines'] = 'Estimate spectral resolution in each slit using the arc lines. '\
                                  'If True, the estimated FWHM will override ``fwhm`` only in '\
                                  'the determination of the wavelength solution (`i.e.`, not in '\
                                  'WaveTilts).'

        # These are the parameters used for reidentification
        defaults['reid_arxiv']=None
        dtypes['reid_arxiv'] = str
        descr['reid_arxiv'] = 'Name of the archival wavelength solution file that will be used ' \
                              'for the wavelength reidentification.  Only used if ``method`` is ' \
                              '\'reidentify\' or \'full_template\'.'

        defaults['nreid_min'] = 1
        dtypes['nreid_min'] = int
        descr['nreid_min'] = 'Minimum number of times that a given candidate reidentified line ' \
                             'must be properly matched with a line in the arxiv to be ' \
                             'considered a good reidentification. If there is a lot of ' \
                             'duplication in the arxiv of the spectra in question (i.e. ' \
                             'multislit) set this to a number like 1-4. For echelle this ' \
                             'depends on the number of solutions in the arxiv.  Set this to 1 ' \
                             'for fixed format echelle spectrographs.  For an echelle with a ' \
                             'tiltable grating, this will depend on the number of solutions in ' \
                             'the arxiv.'

        defaults['wvrng_arxiv'] = None
        dtypes['wvrng_arxiv'] = list
        descr['wvrng_arxiv'] = 'Cut the arxiv template down to this specified wavelength range [min,max]'

        defaults['nsnippet'] = 2
        dtypes['nsnippet'] = int
        descr['nsnippet'] = 'Number of spectra to chop the arc spectrum into when ``method`` is ' \
                            '\'full_template\''

        defaults['cc_thresh'] = 0.70
        dtypes['cc_thresh'] = [float, list, np.ndarray]
        descr['cc_thresh'] = 'Threshold for the *global* cross-correlation coefficient between ' \
                             'an input spectrum and member of the archive required to attempt ' \
                             'reidentification.  Spectra from the archive with a lower ' \
                             'cross-correlation are not used for reidentification. This can be ' \
                             'a single number or a list/array providing the value for each slit.'

        defaults['cc_local_thresh'] = 0.70
        dtypes['cc_local_thresh'] = float
        descr['cc_local_thresh'] = 'Threshold for the *local* cross-correlation coefficient, ' \
                                   'evaluated at each reidentified line,  between an input ' \
                                   'spectrum and the shifted and stretched archive spectrum ' \
                                   'above which a line must be to be considered a good line ' \
                                   'for reidentification. The local cross-correlation is ' \
                                   'evaluated at each candidate reidentified line (using a ' \
                                   'window of nlocal_cc), and is then used to score the the ' \
                                   'reidentified lines to arrive at the final set of good ' \
                                   'reidentifications.'

        defaults['nlocal_cc'] = 11
        dtypes['nlocal_cc'] = int
        descr['nlocal_cc'] = 'Size of pixel window used for local cross-correlation ' \
                             'computation for each arc line. If not an odd number one will ' \
                             'be added to it to make it odd.'

        # These are the parameters used for the iterative fitting of the arc lines
        defaults['rms_threshold'] = 0.15
        dtypes['rms_threshold'] = [float, list, np.ndarray]
        descr['rms_threshold'] = 'Minimum RMS for keeping a slit/order solution. This can be a ' \
                                 'single number or a list/array providing the value for each slit. ' \
                                 'Only used if ``method`` is either \'holy-grail\' or \'reidentify\''

        defaults['match_toler'] = 2.0
        dtypes['match_toler'] = float
        descr['match_toler'] = 'Matching tolerance in pixels when searching for new lines. This ' \
                               'is the difference in pixels between the wavlength assigned to ' \
                               'an arc line by an iteration of the wavelength solution to the ' \
                               'wavelength in the line list.  This parameter is also used as ' \
                               'the matching tolerance in pixels for a line reidentification.  ' \
                               'A good line match must match within this tolerance to the ' \
                               'shifted and stretched archive spectrum, and the archive ' \
                               'wavelength solution at this match must be within match_toler ' \
                               'dispersion elements from the line in line list.'

        defaults['func'] = 'legendre'
        dtypes['func'] = str
        descr['func'] = 'Function used for wavelength solution fits'

        defaults['n_first'] = 2
        dtypes['n_first'] = int
        descr['n_first'] = 'Order of first guess fit to the wavelength solution.'

        defaults['sigrej_first'] = 2.0
        dtypes['sigrej_first'] = float
        descr['sigrej_first'] = 'Number of sigma for rejection for the first guess to the ' \
                                'wavelength solution.'

        defaults['n_final'] = 4
        dtypes['n_final'] = [int, float, list, np.ndarray]
        descr['n_final'] = 'Order of final fit to the wavelength solution (there are n_final+1 ' \
                           'parameters in the fit). This can be a single number or a ' \
                           'list/array providing the value for each slit'

        defaults['sigrej_final'] = 3.0
        dtypes['sigrej_final'] = float
        descr['sigrej_final'] = 'Number of sigma for rejection for the final guess to the ' \
                                'wavelength solution.'

        defaults['numsearch'] = 20
        dtypes['numsearch'] = int
        descr['numsearch'] = 'Number of brightest arc lines to search for in preliminary ' \
                             'identification'

        defaults['nfitpix'] = 5
        dtypes['nfitpix'] = int
        descr['nfitpix'] = 'Number of pixels to fit when deriving the centroid of the arc ' \
                           'lines (an odd number is best)'

        # TODO: What should the default be?  None or 'heliocentric'?
        defaults['refframe'] = 'heliocentric'
        options['refframe'] = WavelengthSolutionPar.valid_reference_frames()
        dtypes['refframe'] = str
        descr['refframe'] = 'Frame of reference for the wavelength calibration.  ' \
                         'Options are: {0}'.format(', '.join(options['refframe']))

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
        k = np.array([*cfg.keys()])
        parkeys = ['reference', 'method', 'echelle', 'ech_nspec_coeff',
                   'ech_norder_coeff', 'ech_sigrej', 'ech_separate_2d', 'lamps', 'sigdetect',
                   'fwhm', 'fwhm_fromlines', 'reid_arxiv', 'nreid_min', 'cc_thresh', 'cc_local_thresh',
                   'nlocal_cc', 'rms_threshold', 'match_toler', 'func', 'n_first','n_final',
                   'sigrej_first', 'sigrej_final', 'numsearch', 'nfitpix',
                   'refframe', 'nsnippet', 'use_instr_flag',
                   'wvrng_arxiv']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for WavelengthSolutionPar.'.format(
                             k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_reference():
        """
        Return the valid wavelength solution methods.
        """
        return ['arc', 'sky', 'pixel']

    @staticmethod
    def valid_methods():
        """
        Return the valid wavelength solution methods.
        """
        return ['holy-grail', 'identify', 'reidentify', 'echelle', 'full_template']

    @staticmethod
    def valid_lamps():
        """
        Return the valid lamp ions
        """
        return ['ArI', 'CdI', 'HgI', 'HeI', 'KrI', 'NeI', 'XeI', 'ZnI', 'ThAr']

    @staticmethod
    def valid_media():
        """
        Return the valid media for the wavelength calibration.
        """
        return ['vacuum', 'air']

    @staticmethod
    def valid_reference_frames():
        """
        Return the valid reference frames for the wavelength calibration
        """
        return ['observed', 'heliocentric', 'barycentric']

    def validate(self):
        pass


class EdgeTracePar(ParSet):
    """
    Parameters used for slit edge tracing.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    prefix = 'ETP'  # Prefix for writing parameters to a header is a class attribute
    def __init__(self, filt_iter=None, sobel_mode=None, edge_thresh=None, sobel_enhance=None, exclude_regions=None,
                 follow_span=None, det_min_spec_length=None, max_shift_abs=None, max_shift_adj=None,
                 max_spat_error=None, match_tol=None, fit_function=None, fit_order=None,
                 fit_maxdev=None, fit_maxiter=None, fit_niter=None, fit_min_spec_length=None,
                 auto_pca=None, left_right_pca=None, pca_min_edges=None, pca_n=None,
                 pca_var_percent=None, pca_function=None, pca_order=None, pca_sigrej=None,
                 pca_maxrej=None, pca_maxiter=None, smash_range=None, edge_detect_clip=None,
                 trace_median_frac=None, trace_thresh=None, fwhm_uniform=None, niter_uniform=None,
                 fwhm_gaussian=None, niter_gaussian=None, det_buffer=None, max_nudge=None,
                 sync_predict=None, sync_center=None, gap_offset=None, sync_to_edge=None,
                 bound_detector=None, minimum_slit_dlength=None, dlength_range=None, 
                 minimum_slit_length=None, minimum_slit_length_sci=None,
                 length_range=None, minimum_slit_gap=None, clip=None, order_match=None,
                 order_offset=None, overlap=None, use_maskdesign=None, maskdesign_maxsep=None,
                 maskdesign_step=None, maskdesign_sigrej=None, pad=None, add_slits=None,
                 add_predict=None, rm_slits=None,
                 maskdesign_filename=None):

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
        defaults['filt_iter'] = 0
        dtypes['filt_iter'] = int
        descr['filt_iter'] = 'Number of median-filtering iterations to perform on sqrt(trace) ' \
                             'image before applying to Sobel filter to detect slit/order edges.'

        defaults['sobel_mode'] = 'nearest'
        options['sobel_mode'] = EdgeTracePar.valid_sobel_modes()
        dtypes['sobel_mode'] = str
        descr['sobel_mode'] = 'Mode for Sobel filtering.  Default is \'nearest\'; note we find' \
                              '\'constant\' works best for DEIMOS.'

        defaults['edge_thresh'] = 20.
        dtypes['edge_thresh'] = [int, float]
        descr['edge_thresh'] = 'Threshold for finding edges in the Sobel-filtered significance image.'

        defaults['sobel_enhance'] = 0
        dtypes['sobel_enhance'] = int
        descr['sobel_enhance'] = 'Enhance the sobel filtering? A value of 0 will not enhance the sobel filtering. ' \
                                 'Any other value > 0 will sum the sobel values. For example, a value of 3 will ' \
                                 'combine the sobel values for the 3 nearest pixels. This is useful when a slit ' \
                                 'edge is poorly defined (e.g. vignetted).'

        defaults['exclude_regions'] = None
        dtypes['exclude_regions'] = [list, str]
        descr['exclude_regions'] = 'User-defined regions to exclude from the slit tracing. To set this parameter, ' \
                                   'the text should be a comma separated list of pixel ranges (in the x direction) ' \
                                   'to be excluded and the detector number. For example, the following string ' \
                                   '1:0:20,1:300:400  would select two regions in det=1 between pixels 0 and 20 ' \
                                   'and between 300 and 400.'

        defaults['follow_span'] = 20
        dtypes['follow_span'] = int
        descr['follow_span'] = 'In the initial connection of spectrally adjacent edge ' \
                               'detections, this sets the number of previous spectral rows ' \
                               'to consider when following slits forward.'

        # TODO: Allow this to be a list so that it can be detector specific?
        defaults['det_min_spec_length'] = 0.33
        dtypes['det_min_spec_length'] = [int, float]
        descr['det_min_spec_length'] = 'The minimum spectral length (as a fraction of the ' \
                                       'detector size) of a trace determined by direct ' \
                                       'measurements of the detector data (as opposed to what ' \
                                       'should be included in any modeling approach; see ' \
                                       'fit_min_spec_length).'

        defaults['max_shift_abs'] = 0.5
        dtypes['max_shift_abs'] = [int, float]
        descr['max_shift_abs'] = 'Maximum spatial shift in pixels between an input edge ' \
                                 'location and the recentroided value.'

        defaults['max_shift_adj'] = 0.15
        dtypes['max_shift_adj'] = [int, float]
        descr['max_shift_adj'] = 'Maximum spatial shift in pixels between the edges in ' \
                                 'adjacent spectral positions.'

#        defaults['max_spat_error'] = 0.2
        dtypes['max_spat_error'] = [int, float]
        descr['max_spat_error'] = 'Maximum error in the spatial position of edges in pixels.'

        defaults['match_tol'] = 3.
        dtypes['match_tol'] = [int, float]
        descr['match_tol'] = 'Same-side slit edges below this separation in pixels are ' \
                             'considered part of the same edge.'

        defaults['fit_function'] = 'legendre'
        options['fit_function'] = EdgeTracePar.valid_functions()
        dtypes['fit_function'] = str
        descr['fit_function'] = 'Function fit to edge measurements.  ' \
                                'Options are: {0}'.format(', '.join(options['fit_function']))

        defaults['fit_order'] = 5
        dtypes['fit_order'] = int
        descr['fit_order'] = 'Order of the function fit to edge measurements.'

        defaults['fit_maxdev'] = 5.0
        dtypes['fit_maxdev'] = [int, float]
        descr['fit_maxdev'] = 'Maximum deviation between the fitted and measured edge position ' \
                              'for rejection in spatial pixels.'

        defaults['fit_maxiter'] = 25
        dtypes['fit_maxiter'] = int
        descr['fit_maxiter'] = 'Maximum number of rejection iterations during edge fitting.'

        defaults['fit_niter'] = 1
        dtypes['fit_niter'] = int
        descr['fit_niter'] = 'Number of iterations of re-measuring and re-fitting the edge ' \
                             'data; see :func:`pypeit.core.trace.fit_trace`.'

        # TODO: Allow this to be a list so that it can be detector specific?
        defaults['fit_min_spec_length'] = 0.6
        dtypes['fit_min_spec_length'] = float
        descr['fit_min_spec_length'] = 'Minimum unmasked spectral length of a traced slit edge ' \
                                       'to use in any modeling procedure (polynomial fitting ' \
                                       'or PCA decomposition).'

        defaults['auto_pca'] = True
        dtypes['auto_pca'] = bool
        descr['auto_pca'] = 'During automated tracing, attempt to construct a PCA decomposition ' \
                            'of the traces. When True, the edge traces resulting from the ' \
                            'initial detection, centroid refinement, and polynomial fitting ' \
                            'must meet a set of criteria for performing the pca; see ' \
                            ':func:`pypeit.edgetrace.EdgeTraceSet.can_pca`.  If False, the ' \
                            '``sync_predict`` parameter *cannot* be set to ``pca``; if it is ' \
                            'not, the value is set to ``nearest`` and a warning is issued when ' \
                            'validating the parameter set.'

        defaults['left_right_pca'] = False
        dtypes['left_right_pca'] = bool
        descr['left_right_pca'] = 'Construct a PCA decomposition for the left and right traces ' \
                                  'separately.  This can be important for cross-dispersed ' \
                                  'echelle spectrographs (e.g., Keck-NIRES)'

        defaults['pca_min_edges'] = 4
        dtypes['pca_min_edges'] = int
        descr['pca_min_edges'] = 'Minimum number of edge traces required to perform a PCA ' \
                                 'decomposition of the trace form.  If left_right_pca is True, ' \
                                 'this minimum applies to the number of left and right traces ' \
                                 'separately.'

        dtypes['pca_n'] = int
        descr['pca_n'] = 'The number of PCA components to keep, which must be less than the ' \
                         'number of detected traces.  If not provided, determined by ' \
                         'calculating the minimum number of components required to explain a ' \
                         'given percentage of variance in the edge data; see `pca_var_percent`.'

        defaults['pca_var_percent'] = 99.8
        dtypes['pca_var_percent'] = [int, float]
        descr['pca_var_percent'] = 'The percentage (i.e., not the fraction) of the variance in ' \
                                   'the edge data accounted for by the PCA used to truncate ' \
                                   'the number of PCA coefficients to keep (see `pca_n`).  ' \
                                   'Ignored if `pca_n` is provided directly.'

        defaults['pca_function'] = 'polynomial'
        dtypes['pca_function'] = str
        options['pca_function'] = EdgeTracePar.valid_functions()
        descr['pca_function'] = 'Type of function fit to the PCA coefficients for each ' \
                                'component.  Options are: {0}'.format(
                                    ', '.join(options['pca_function']))

        defaults['pca_order'] = 2
        dtypes['pca_order'] = int
        descr['pca_order'] = 'Order of the function fit to the PCA coefficients.'

        defaults['pca_sigrej'] = [2., 2.]
        dtypes['pca_sigrej'] = [int, float, list]
        descr['pca_sigrej'] = 'Sigma rejection threshold for fitting PCA components. Individual ' \
                              'numbers are used for both lower and upper rejection. A list of ' \
                              'two numbers sets these explicitly (e.g., [2., 3.]).'

        defaults['pca_maxrej'] = 1
        dtypes['pca_maxrej'] = int
        descr['pca_maxrej'] = 'Maximum number of PCA coefficients rejected during a given fit ' \
                              'iteration.'

        defaults['pca_maxiter'] = 25
        dtypes['pca_maxiter'] = int
        descr['pca_maxiter'] = 'Maximum number of rejection iterations when fitting the PCA ' \
                               'coefficients.'

        defaults['smash_range'] = [0., 1.]
        dtypes['smash_range'] = list
        descr['smash_range'] = 'Range of the slit in the spectral direction (in fractional ' \
                               'units) to smash when searching for slit edges.  If the ' \
                               'spectrum covers only a portion of the image, use that range.'

        # TODO: Does this still need to be different from `edge_thresh`?
        dtypes['edge_detect_clip'] = [int, float]
        descr['edge_detect_clip'] = 'Sigma clipping level for peaks detected in the collapsed, ' \
                                    'Sobel-filtered significance image.'

        dtypes['trace_median_frac'] = [int, float]
        descr['trace_median_frac'] = 'After detection of peaks in the rectified Sobel-filtered ' \
                                     'image and before refitting the edge traces, the rectified ' \
                                     'image is median filtered with a kernel width of ' \
                                     '`trace_median_frac*nspec` along the spectral dimension.'

        dtypes['trace_thresh'] = [int, float]
        descr['trace_thresh'] = 'After rectification and median filtering of the Sobel-filtered ' \
                                'image (see `trace_median_frac`), values in the median-filtered ' \
                                'image *below* this threshold are masked in the refitting of ' \
                                'the edge trace data.  If None, no masking applied.'

        defaults['fwhm_uniform'] = 3.0
        dtypes['fwhm_uniform'] = [int, float]
        descr['fwhm_uniform'] = 'The `fwhm` parameter to use when using uniform weighting in ' \
                                ':func:`pypeit.core.trace.fit_trace` when refining the PCA ' \
                                'predictions of edges.  See description of ' \
                                ':func:`pypeit.core.trace.peak_trace`.'

        defaults['niter_uniform'] = 9
        dtypes['niter_uniform'] = int
        descr['niter_uniform'] = 'The number of iterations of ' \
                                 ':func:`pypeit.core.trace.fit_trace` to use when using ' \
                                 'uniform weighting.'

        defaults['fwhm_gaussian'] = 3.0
        dtypes['fwhm_gaussian'] = [int, float]
        descr['fwhm_gaussian'] = 'The `fwhm` parameter to use when using Gaussian weighting in ' \
                                 ':func:`pypeit.core.trace.fit_trace` when refining the PCA ' \
                                 'predictions of edges.  See description ' \
                                 ':func:`pypeit.core.trace.peak_trace`.'

        defaults['niter_gaussian'] = 6
        dtypes['niter_gaussian'] = int
        descr['niter_gaussian'] = 'The number of iterations of ' \
                                  ':func:`pypeit.core.trace.fit_trace` to use when using ' \
                                  'Gaussian weighting.'

        defaults['det_buffer'] = 5
        dtypes['det_buffer'] = int
        descr['det_buffer'] = 'The minimum separation between the detector edges and a slit ' \
                              'edge for any added edge traces.  Must be positive.'

#        defaults['max_nudge'] = 100
        dtypes['max_nudge'] = [int, float]
        descr['max_nudge'] = 'If parts of any (predicted) trace fall off the detector edge, ' \
                             'allow them to be nudged away from the detector edge up to and ' \
                             'including this maximum number of pixels.  If None, no limit is ' \
                             'set; otherwise should be 0 or larger.'

        defaults['sync_predict'] = 'pca'
        options['sync_predict'] = EdgeTracePar.valid_predict_modes()
        dtypes['sync_predict'] = str
        descr['sync_predict'] = 'Mode to use when predicting the form of the trace to insert.  ' \
                                'Use `pca` to use the PCA decomposition, `nearest` to ' \
                                'reproduce the shape of the nearest trace, or `auto` to let PypeIt ' \
                                'decide which mode to use between `pca` and `nearest`. In general, ' \
                                'it will first try `pca`, and if that is not possible, it will use `nearest`.'

        defaults['sync_center'] = 'median'
        options['sync_center'] = EdgeTracePar.valid_center_modes()
        dtypes['sync_center'] = str
        descr['sync_center'] = 'Mode to use for determining the location of traces to insert.  ' \
                               'Use `median` to use the median of the matched left and right ' \
                               'edge pairs, `nearest` to use the length of the nearest slit, ' \
                               'or `gap` to offset by a fixed gap width from the next slit edge.'

        defaults['gap_offset'] = 5.
        dtypes['gap_offset'] = [int, float]
        descr['gap_offset'] = 'Offset (pixels) used for the slit edge gap width when inserting ' \
                              'slit edges (see `sync_center`) or when nudging predicted slit ' \
                              'edges to avoid slit overlaps.  This should be larger than ' \
                              '`minimum_slit_gap` when converted to arcseconds.'

        defaults['sync_to_edge'] = True
        dtypes['sync_to_edge'] = bool
        descr['sync_to_edge'] = 'If adding a first left edge or a last right edge, ignore ' \
                                '`center_mode` for these edges and place them at the edge of ' \
                                'the detector (with the relevant shape).'

        defaults['bound_detector'] = False
        dtypes['bound_detector'] = bool
        descr['bound_detector'] = 'When the code is ready to synchronize the left/right trace ' \
                                  'edges, the traces should have been constructed, vetted, and ' \
                                  'cleaned. This can sometimes lead to *no* valid traces. This ' \
                                  'parameter dictates what to do next. If ``bound_detector`` is ' \
                                  'True, the code will artificially add left and right edges ' \
                                  'that bound the detector; if False, the code identifies the ' \
                                  'slit-edge tracing as being unsuccessful, warns the user, and ' \
                                  'ends gracefully. Note that setting ``bound_detector`` to ' \
                                  'True is needed for some long-slit data where the slit ' \
                                  'edges are, in fact, beyond the edges of the detector.'

        dtypes['minimum_slit_dlength'] = [int, float]
        descr['minimum_slit_dlength'] = 'Minimum *change* in the slit length (arcsec) as a ' \
                                        'function of wavelength in arcsec.  This is mostly ' \
                                        'meant to catch cases when the polynomial fit to the ' \
                                        'detected edges becomes ill-conditioned (e.g., when ' \
                                        'the slits run off the edge of the detector) and leads ' \
                                        'to wild traces.  If reducing the order of the ' \
                                        'polynomial (``fit_order``) does not help, try using ' \
                                        'this to remove poorly constrained slits.'

        dtypes['dlength_range'] = [int, float]
        descr['dlength_range'] = 'Similar to ``minimum_slit_dlength``, but constrains the ' \
                                 '*fractional* change in the slit length as a function of ' \
                                 'wavelength.  For example, a value of 0.2 means that slit ' \
                                 'length should not vary more than 20%' \
                                 'as a function of wavelength.  '

#        defaults['minimum_slit_length'] = 6.
        dtypes['minimum_slit_length'] = [int, float]
        descr['minimum_slit_length'] = 'Minimum slit length in arcsec.  Slit lengths are ' \
                                       'determined by the median difference between the left ' \
                                       'and right edge locations for the unmasked trace ' \
                                       'locations.  This is used to identify traces that are ' \
                                       '*erroneously* matched together to form slits.  Short ' \
                                       'slits are expected to be ignored or removed (see ' \
                                       ' ``clip``).  If None, no minimum slit length applied.'

        dtypes['minimum_slit_length_sci'] = [int, float]
        descr['minimum_slit_length_sci'] = 'Minimum slit length in arcsec for a science slit.  ' \
                                       'Slit lengths are determined by the median difference ' \
                                       'between the left and right edge locations for the ' \
                                       'unmasked trace locations.  Used in combination with ' \
                                       '``minimum_slit_length``, this parameter is used to ' \
                                       'identify box or alignment slits; i.e., those slits ' \
                                       'that are shorter than ``minimum_slit_length_sci`` but ' \
                                       'larger than ``minimum_slit_length`` are box/alignment ' \
                                       'slits.  Box slits are *never* removed (see ``clip``), ' \
                                       'but no spectra are extracted from them.  If None, no ' \
                                       'minimum science slit length is applied.'

#        defaults['length_range'] = 0.3
        dtypes['length_range'] = [int, float]
        descr['length_range'] = 'Allowed range in slit length compared to the median slit ' \
                                'length.  For example, a value of 0.3 means that slit lengths ' \
                                'should not vary more than 30%.  Relatively shorter or longer ' \
                                'slits are masked or clipped.  Most useful for echelle or ' \
                                'multi-slit data where the slits should have similar or ' \
                                'identical lengths.'

        # TODO: Define this in pixels instead of arcsec?
        dtypes['minimum_slit_gap'] = [int, float]
        descr['minimum_slit_gap'] = 'Minimum slit gap in arcsec.  Gaps between slits are ' \
                                    'determined by the median difference between the right ' \
                                    'and left edge locations of adjacent slits.  Slits with ' \
                                    'small gaps are merged by removing the intervening traces.' \
                                    'If None, no minimum slit gap is applied.  This should be ' \
                                    'smaller than `gap_offset` when converted to pixels.'

        defaults['clip'] = True
        dtypes['clip'] = bool
        descr['clip'] = 'Remove traces flagged as bad, instead of only masking them.  This ' \
                        'is currently only used by ' \
                        ':func:`~pypeit.edgetrace.EdgeTraceSet.centroid_refine`.'

        dtypes['order_match'] = [int, float]
        descr['order_match'] = 'For echelle spectrographs, this is the tolerance allowed for ' \
                               'matching identified "slits" to echelle orders. Must be in ' \
                               'the fraction of the detector spatial scale (i.e., a value of ' \
                               '0.05 means that the order locations must be within 5% of the ' \
                               'expected value).  If None, no limit is used.'

        dtypes['order_offset'] = [int, float]
        descr['order_offset'] = 'Offset to introduce to the expected order positions to improve ' \
                                'the match for this specific data. This is an additive offset ' \
                                'to the measured slit positions; i.e., this should minimize the ' \
                                'difference between the expected order positions and ' \
                                '``self.slit_spatial_center() + offset``. Must be in the ' \
                                'fraction of the detector spatial scale. If None, no offset ' \
                                'is applied.'

        defaults['overlap'] = False
        dtypes['overlap'] = bool
        descr['overlap'] = 'Assume slits identified as abnormally short are actually due to ' \
                           'overlaps between adjacent slits/orders.  If set to True, you *must* ' \
                           'have also used ``length_range`` to identify left-right edge pairs ' \
                           'that have an abnormally short separation.  For those short slits, ' \
                           'the code attempts to convert the short slits into slit gaps.  This ' \
                           'is particularly useful for blue orders in Keck-HIRES data.'

        defaults['use_maskdesign'] = False
        dtypes['use_maskdesign'] = bool
        descr['use_maskdesign'] = 'Use slit-mask designs to identify slits.'

        defaults['maskdesign_filename'] = None
        dtypes['maskdesign_filename'] = [str, list]
        descr['maskdesign_filename'] = 'Mask design info contained in this file or files (comma separated)'

        defaults['maskdesign_maxsep'] = 50
        dtypes['maskdesign_maxsep'] = [int, float]
        descr['maskdesign_maxsep'] = 'Maximum allowed offset in pixels between the slit edges ' \
                                     'defined by the slit-mask design and the traced edges.'

        defaults['maskdesign_step'] = 1
        dtypes['maskdesign_step'] = [int, float]
        descr['maskdesign_step'] = 'Step in pixels used to generate a list of possible offsets ' \
                                   '(within +/- `maskdesign_maxsep`) between the slit edges defined ' \
                                   'by the mask design and the traced edges.'

        defaults['maskdesign_sigrej'] = 3
        dtypes['maskdesign_sigrej'] = [int, float]
        descr['maskdesign_sigrej'] = 'Number of sigma for sigma-clipping rejection during slit-mask ' \
                                     'design matching.'

#        # Force trim to be a tuple
#        if pars['trim'] is not None and not isinstance(pars['trim'], tuple):
#            try:
#                pars['trim'] = tuple(pars['trim'])
#            except:
#                raise TypeError('Could not convert provided trim to a tuple.')
#        defaults['trim'] = (0,0)
#        dtypes['trim'] = tuple
#        descr['trim'] = 'How much to trim off each edge of each slit.  Each number should be 0 ' \
#                        'or positive'

        # TODO: Describe better where and how this is used.  It's not
        # actually used in the construction of the nominal slit edges,
        # but only in subsequent use of the slits (e.g., flat-fielding)
        defaults['pad'] = 0
        dtypes['pad'] = int
        descr['pad'] = 'Integer number of pixels to consider beyond the slit edges when ' \
                       'selecting pixels that are \'on\' the slit.'

#        defaults['single'] = []
#        dtypes['single'] = list
#        descr['single'] = 'Add a single, user-defined slit based on its location on each ' \
#                          'detector.  Syntax is a list of values, 2 per detector, that define ' \
#                          'the slit according to column values.  The second value (for the ' \
#                          'right edge) must be greater than 0 to be applied.  LRISr example: ' \
#                          'setting single = -1, -1, 7, 295 means the code will skip the ' \
#                          'user-definition for the first detector but adds one for the second. ' \
#                          ' None means no user-level slits defined.'

        dtypes['add_slits'] = [str, list]
        descr['add_slits'] = 'Add one or more user-defined slits.  The syntax to define a ' \
                             'slit to add is: \'det:spec:spat_left:spat_right\' where ' \
                             'det=detector, spec=spectral pixel, spat_left=spatial pixel of ' \
                             'left slit boundary, and spat_righ=spatial pixel of right slit ' \
                             'boundary.  For example, \'2:2000:2121:2322,3:2000:1201:1500\' ' \
                             'will add a slit to detector 2 passing through spec=2000 ' \
                             'extending spatially from 2121 to 2322 and another on detector 3 ' \
                             'at spec=2000 extending from 1201 to 1500.'

        defaults['add_predict'] = 'nearest'
        dtypes['add_predict'] = str
        descr['add_predict'] = 'Sets the method used to predict the shape of the left and right ' \
                               'traces for a user-defined slit inserted.  Options are (1) ' \
                               '``straight`` inserts traces with a constant spatial pixels ' \
                               'position, (2) ``nearest`` inserts traces with a form identical ' \
                               'to the automatically identified trace at the nearest spatial ' \
                               'position to the inserted slit, or (3) ``pca`` uses the PCA ' \
                               'decomposition to predict the shape of the traces.'

        dtypes['rm_slits'] = [str, list]
        descr['rm_slits'] = 'Remove one or more user-specified slits.  The syntax used to ' \
                            'define a slit to remove is: \'det:spec:spat\' where det=detector, ' \
                            'spec=spectral pixel, spat=spatial pixel.  For example, ' \
                            '\'2:2000:2121,3:2000:1500\' will remove the slit on detector 2 ' \
                            'that contains pixel (spat,spec)=(2000,2121) and on detector 3 ' \
                            'that contains pixel (2000,2121).'

        # Instantiate the parameter set
        super(EdgeTracePar, self).__init__(list(pars.keys()), values=list(pars.values()),
                                           defaults=list(defaults.values()),
                                           options=list(options.values()),
                                           dtypes=list(dtypes.values()),
                                           descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        # TODO Please provide docs
        k = np.array([*cfg.keys()])
        parkeys = ['filt_iter', 'sobel_mode', 'edge_thresh', 'sobel_enhance', 'exclude_regions',
                   'follow_span', 'det_min_spec_length',
                   'max_shift_abs', 'max_shift_adj', 'max_spat_error', 'match_tol', 'fit_function',
                   'fit_order', 'fit_maxdev', 'fit_maxiter', 'fit_niter', 'fit_min_spec_length',
                   'auto_pca', 'left_right_pca', 'pca_min_edges', 'pca_n', 'pca_var_percent',
                   'pca_function', 'pca_order', 'pca_sigrej', 'pca_maxrej', 'pca_maxiter',
                   'smash_range', 'edge_detect_clip', 'trace_median_frac', 'trace_thresh',
                   'fwhm_uniform', 'niter_uniform', 'fwhm_gaussian', 'niter_gaussian',
                   'det_buffer', 'max_nudge', 'sync_predict', 'sync_center', 'gap_offset',
                   'sync_to_edge', 'bound_detector', 'minimum_slit_dlength', 'dlength_range', 
                   'minimum_slit_length', 'minimum_slit_length_sci', 'length_range',
                   'minimum_slit_gap', 'clip', 'order_match', 'order_offset', 'overlap',
                   'use_maskdesign', 'maskdesign_maxsep', 'maskdesign_step', 
                   'maskdesign_sigrej', 'maskdesign_filename',
                   'pad', 'add_slits', 'add_predict', 'rm_slits']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for EdgeTracePar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_functions():
        """
        Return the list of valid functions to use for slit tracing.
        """
        return ['polynomial', 'legendre', 'chebyshev']

    @staticmethod
    def valid_sobel_modes():
        """Return the valid sobel modes."""
        return ['nearest', 'constant']

    @staticmethod
    def valid_predict_modes():
        """Return the valid trace prediction modes."""
        return ['pca', 'nearest', 'auto']

    @staticmethod
    def valid_center_modes():
        """Return the valid center prediction modes."""
        return ['median', 'nearest', 'gap']

    def validate(self):
        """Validate the parameter set."""
        if not self['auto_pca'] and self['sync_predict'] == 'pca':
            warnings.warn('sync_predict cannot be pca if auto_pca is False.  Setting to nearest.')
            self['sync_predict'] = 'nearest'


class WaveTiltsPar(ParSet):
    """
    The parameter set used to hold arguments for tracing the
    monochromatic tilt along the slit.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.

    .. todo::
        Changed to reflect wavetilts.py settings.  Was `yorder`
        previously `disporder`?  If so, I think I prefer the generality
        of `disporder`...
    """
    def __init__(self, idsonly=None, tracethresh=None, sig_neigh=None, nfwhm_neigh=None,
                 maxdev_tracefit=None, sigrej_trace=None, spat_order=None, spec_order=None,
                 func2d=None, maxdev2d=None, sigrej2d=None, rm_continuum=None, cont_rej=None,
                 minmax_extrap=None):

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
                           'tilts (CURRENTLY NOT USED!)'

        defaults['tracethresh'] = 20.
        dtypes['tracethresh'] = [int, float, list, np.ndarray]
        descr['tracethresh'] = 'Significance threshold for arcs to be used in tracing wavelength tilts. ' \
                               'This can be a single number or a list/array providing the value for each slit/order.'


        defaults['sig_neigh'] = 10.
        dtypes['sig_neigh'] = [int, float]
        descr['sig_neigh'] = 'Significance threshold for arcs to be used in line identification for the purpose of identifying neighboring lines. ' \
                             'The tracethresh parameter above determines the significance threshold of lines that will be traced, but these lines ' \
                             ' must be at least nfwhm_neigh fwhm away from neighboring lines. This parameter determines the significance above which ' \
                             ' a line must be to be considered a possible colliding neighbor. A low value of sig_neigh will result in an overall ' \
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
        dtypes['spat_order'] = [int, float, list, np.ndarray]
        descr['spat_order'] = 'Order of the legendre polynomial to be fit to the the tilt of an arc line. This parameter determines ' \
                              'both the orer of the *individual* arc line tilts, as well as the order of the spatial direction of the ' \
                              '2d legendre polynomial (spatial, spectral) that is fit to obtain a global solution for the tilts across the ' \
                              'slit/order. This can be a single number or a list/array providing the value for each slit'

        defaults['spec_order'] = 4
        dtypes['spec_order'] = [int, float, list, np.ndarray]
        descr['spec_order'] = 'Order of the spectral direction of the 2d legendre polynomial (spatial, spectral) that is ' \
                              'fit to obtain a global solution for the tilts across the slit/order. ' \
                              'This can be a single number or a list/array providing the value for each slit'


        defaults['minmax_extrap'] = [150., 1000.]
        dtypes['minmax_extrap'] = [list, np.ndarray]
        descr['minmax_extrap'] = 'Sets how far below the last measured tilt line is extrapolated in tracewave.fit_tilts()'

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

        defaults['rm_continuum'] = False
        dtypes['rm_continuum'] = bool
        descr['rm_continuum'] = 'Before tracing the line center at each spatial position, ' \
                                'remove any low-order continuum in the 2D spectra.'

        # TODO: Replace these with relevant parameters from
        # arc.iter_continuum
#        defaults['cont_function'] = 'legendre'
#        dtypes['cont_function'] = str
#        descr['cont_function'] = 'Function type used to fit the continuum to be removed.'
#
#        defaults['cont_order'] = 3
#        dtypes['cont_order'] = int
#        descr['cont_order'] = 'Order of the function used to fit the continuum to be removed.'

        defaults['cont_rej'] = [3, 1.5]
        dtypes['cont_rej'] = [int, float, list, np.ndarray]
        descr['cont_rej'] = 'The sigma threshold for rejection.  Can be a single number or two ' \
                            'numbers that give the low and high sigma rejection, respectively.'

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
        k = np.array([*cfg.keys()])
        parkeys = ['idsonly', 'tracethresh', 'sig_neigh', 'maxdev_tracefit', 'sigrej_trace',
                   'nfwhm_neigh', 'spat_order', 'spec_order', 'func2d', 'maxdev2d', 'sigrej2d',
                   'rm_continuum', 'cont_rej', 'minmax_extrap'] #'cont_function', 'cont_order',

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for WaveTiltsPar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)


    def validate(self):
        if hasattr(self.data['cont_rej'], '__len__'):
            if len(self.data['cont_rej']) != 2:
                raise ValueError('Continuum rejection threshold must be a single number or a '
                                 'two-element list/array.')


class ReducePar(ParSet):
    """
    The parameter set used to hold arguments for sky subtraction, object
    finding and extraction in the Reduce class

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """

    def __init__(self, findobj=None, skysub=None, extraction=None,
                 cube=None, trim_edge=None, slitmask=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k, values[k]) for k in args[1:]])  # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['findobj'] = FindObjPar()
        dtypes['findobj'] = [ ParSet, dict ]
        descr['findobj'] = 'Parameters for the find object and tracing algorithms'

        defaults['skysub'] = SkySubPar()
        dtypes['skysub'] = [ ParSet, dict ]
        descr['skysub'] = 'Parameters for sky subtraction algorithms'

        defaults['extraction'] = ExtractionPar()
        dtypes['extraction'] = [ ParSet, dict ]
        descr['extraction'] = 'Parameters for extraction algorithms'

        defaults['slitmask'] = SlitMaskPar()
        dtypes['slitmask'] = [ ParSet, dict ]
        descr['slitmask'] = 'Parameters for slitmask'

        defaults['cube'] = CubePar()
        dtypes['cube'] = [ ParSet, dict ]
        descr['cube'] = 'Parameters for cube generation algorithms'

        defaults['trim_edge'] = [3, 3]
        dtypes['trim_edge'] = list
        descr['trim_edge'] = 'Trim the slit by this number of pixels left/right when performing sky subtraction'

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
        k = np.array([*cfg.keys()])

        allkeys = ['findobj', 'skysub', 'extraction', 'cube', 'trim_edge', 'slitmask']
        badkeys = np.array([pk not in allkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for ReducePar.'.format(k[badkeys]))

        kwargs = {}
        # Keywords that are ParSets
        pk = 'findobj'
        kwargs[pk] = FindObjPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'skysub'
        kwargs[pk] = SkySubPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'extraction'
        kwargs[pk] = ExtractionPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'cube'
        kwargs[pk] = CubePar.from_dict(cfg[pk]) if pk in k else None
        pk = 'slitmask'
        kwargs[pk] = SlitMaskPar.from_dict(cfg[pk]) if pk in k else None

        return cls(**kwargs)

    def validate(self):
        pass


class FindObjPar(ParSet):
    """
    The parameter set used to hold arguments for functionality relevant
    to finding and tracing objects.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """

    def __init__(self, trace_npoly=None, snr_thresh=None, find_trim_edge=None,
                 find_maxdev=None, find_extrap_npoly=None, maxnumber_sci=None, maxnumber_std=None,
                 find_fwhm=None, ech_find_max_snr=None, ech_find_min_snr=None,
                 ech_find_nabove_min_snr=None, skip_second_find=None, skip_final_global=None,
                 skip_skysub=None, find_negative=None, find_min_max=None, std_spec1d=None):
        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k, values[k]) for k in args[1:]])  # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['trace_npoly'] = 5
        dtypes['trace_npoly'] = int
        descr['trace_npoly'] = 'Order of legendre polynomial fits to object traces.'

        defaults['maxnumber_sci'] = 10
        dtypes['maxnumber_sci'] = int
        descr['maxnumber_sci'] = 'Maximum number of objects to extract in a science frame.  Use ' \
                             'None for no limit. This parameter can be useful in situations where systematics lead to ' \
                             'spurious extra objects. Setting this parameter means they will be trimmed. ' \
                             'For mulitslit maxnumber applies per slit, for echelle observations this ' \
                             'applies per order. Note that objects on a slit/order impact the sky-modeling and so ' \
                             'maxnumber should never be lower than the true number of detectable objects on your slit. ' \
                             'For image differenced observations with positive and negative object traces, maxnumber applies ' \
                             'to the number of positive (or negative) traces individually. In other words, if you had two positive objects and ' \
                             'one negative object, then you would set maxnumber to be equal to two (not three). Note that if manually ' \
                             'extracted apertures are explicitly requested, they do not count against this maxnumber. If more than ' \
                             'maxnumber objects are detected, then highest S/N ratio objects will be the ones that are kept. ' \
                             'For multislit observations the choice here depends on the slit length. For echelle observations ' \
                             'with short slits we set the default to be 1'

        defaults['maxnumber_std'] = 5
        dtypes['maxnumber_std'] = int
        descr['maxnumber_std'] = 'Maximum number of objects to extract in a standard star frame.  Same functionality as ' \
                                 'maxnumber_sci documented above. For multislit observations the default here is 5, for echelle ' \
                                 'observations the default is 1'

        defaults['snr_thresh'] = 10.0
        dtypes['snr_thresh'] = [int, float]
        descr['snr_thresh'] = 'S/N threshold for object finding in wavelength direction smashed image.'

        defaults['find_trim_edge'] = [5,5]
        dtypes['find_trim_edge'] = list
        descr['find_trim_edge'] = 'Trim the slit by this number of pixels left/right before finding objects'

        defaults['find_extrap_npoly'] = 3
        dtypes['find_extrap_npoly'] = int
        descr['find_extrap_npoly'] = 'Polynomial order used for trace extrapolation'

        defaults['find_maxdev'] = 2.0
        dtypes['find_maxdev'] = [int, float]
        descr['find_maxdev'] = 'Maximum deviation of pixels from polynomial fit to trace used to reject bad pixels in trace fitting.'

        defaults['find_fwhm'] = 5.0
        dtypes['find_fwhm'] = [int, float]
        descr['find_fwhm'] = 'Indicates roughly the fwhm of objects in pixels for object finding'

        defaults['ech_find_max_snr'] = 1.0
        dtypes['ech_find_max_snr'] = [int, float]
        descr['ech_find_max_snr'] = 'Criteria for keeping echelle objects. They must either have a maximum S/N across all the orders greater than this value ' \
                                    ' or satisfy the min_snr criteria described by the min_snr parameters. If maxnumber is set (see above) then these criteria ' \
                                    'will be applied but only the maxnumber highest (median) S/N ratio objects will be kept. '

        defaults['ech_find_min_snr'] = 0.3
        dtypes['ech_find_min_snr'] = [int, float]
        descr['ech_find_min_snr'] = 'Criteria for keeping echelle objects. They must either have a maximum S/N across all the orders greater than ech_find_max_snr,  value ' \
                                    ' or they must have S/N > ech_find_min_snr on >= ech_find_nabove_min_snr orders. If maxnumber is set (see above) then these criteria ' \
                                    'will be applied but only the maxnumber highest (median) S/N ratio objects will be kept. '

        defaults['ech_find_nabove_min_snr'] = 2
        dtypes['ech_find_nabove_min_snr'] = int
        descr['ech_find_nabove_min_snr'] = 'Criteria for keeping echelle objects. They must either have a maximum S/N across ' \
                                           'all the orders greater than ech_find_max_snr,  value ' \
                                           ' or they must have S/N > ech_find_min_snr on >= ech_find_nabove_min_snr orders. ' \
                                           'If maxnumber is set (see above) then these criteria ' \
                                           'will be applied but only the maxnumber highest (median) S/N ratio objects will be kept.'


        defaults['skip_second_find'] = False
        dtypes['skip_second_find'] = bool
        descr['skip_second_find'] = 'Only perform one round of object finding (mainly for quick_look)'

        defaults['skip_final_global'] = False
        dtypes['skip_final_global'] = bool
        descr['skip_final_global'] = 'If True, do not update initial sky to get global sky using updated noise model. This ' \
                                     'should be True for quicklook to save time. This should also be True for near-IR ' \
                                     'reductions which perform difference imaging, since there we fit sky-residuals rather ' \
                                     'than the sky itself, so there is no noise model to update. '

        defaults['skip_skysub'] = False
        dtypes['skip_skysub'] = bool
        descr['skip_skysub'] = 'If True, do not sky subtract when performing object finding. This should be set to ' \
                               'True for example when running on data that is already sky-subtracted. ' \
                               'Note that for near-IR difference imaging one still wants to remove sky-residuals via ' \
                               'sky-subtraction, and so this is typically set to False'



        defaults['find_negative'] = None
        dtypes['find_negative'] = bool
        descr['find_negative'] = 'Identify negative objects in object finding for spectra that are differenced. This is used to manually ' \
                                 'override the default behavior in PypeIt for object finding by setting this parameter to something other than None ' \
                                 'The default behavior is that PypeIt will search for negative object traces if background frames ' \
                                 'are present in the PypeIt file that are classified as "science" ' \
                                 '(i.e. via pypeit_setup -b, and setting bkg_id in the PypeIt file). If background frames are present ' \
                                 'that are classified as "sky", then PypeIt will NOT search for negative object traces. If one wishes ' \
                                 'to explicitly override this default behavior, set this parameter to True to find negative objects or False to ignore ' \
                                 'them.'

        defaults['find_min_max'] = None
        dtypes['find_min_max'] = list
        descr['find_min_max'] = 'It defines the minimum and maximum of your object in the spectral direction on the ' \
                                'detector. It only used for object finding. This parameter is helpful if your object only ' \
                                'has emission lines or at high redshift and the trace only shows in part of the detector.'

        defaults['std_spec1d'] = None
        dtypes['std_spec1d'] = str
        descr['std_spec1d'] = 'A PypeIt spec1d file of a previously reduced standard star.  The ' \
                              'trace of the standard star spectrum is used as a crutch for ' \
                              'tracing the object spectra, when a direct trace is not possible ' \
                              '(i.e., faint sources).  If provided, this overrides use of any ' \
                              'standards included in your pypeit file; the standard exposures ' \
                              'will still be reduced.'

        # Instantiate the parameter set
        super(FindObjPar, self).__init__(list(pars.keys()),
                                        values=list(pars.values()),
                                        defaults=list(defaults.values()),
                                        options=list(options.values()),
                                        dtypes=list(dtypes.values()),
                                        descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])

        # Basic keywords
        parkeys = ['trace_npoly', 'snr_thresh', 'find_trim_edge',
                   'find_extrap_npoly', 'maxnumber_sci', 'maxnumber_std',
                   'find_maxdev', 'find_fwhm', 'ech_find_max_snr',
                   'ech_find_min_snr', 'ech_find_nabove_min_snr', 'skip_second_find',
                   'skip_final_global', 'skip_skysub', 'find_negative', 'find_min_max',
                   'std_spec1d']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for FindObjPar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        if self.data['std_spec1d'] is not None \
                and not Path(self.data['std_spec1d']).resolve().exists():
            msgs.error(f'{self.data["std_spec1d"]} does not exist!')


class SkySubPar(ParSet):
    """
    The parameter set used to hold arguments for functionality relevant
    to sky subtraction.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """

    def __init__(self, bspline_spacing=None, sky_sigrej=None, global_sky_std=None, no_poly=None,
                 user_regions=None, joint_fit=None, mask_by_boxcar=None,
                 no_local_sky=None, max_mask_frac=None, local_maskwidth=None):
        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k, values[k]) for k in args[1:]])  # "1:" to skip 'self'

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



        defaults['global_sky_std'] = True
        dtypes['global_sky_std'] = bool
        descr['global_sky_std'] = 'Global sky subtraction will be performed on standard stars. This should be turned ' \
                                  'off for example for near-IR reductions with narrow slits, since bright standards can ' \
                                  'fill the slit causing global sky-subtraction to fail. In these situations we go ' \
                                  'straight to local sky-subtraction since it is designed to deal with such situations'

        defaults['no_poly'] = False
        dtypes['no_poly'] = bool
        descr['no_poly'] = 'Turn off polynomial basis (Legendre) in global sky subtraction'

        defaults['no_local_sky'] = False
        dtypes['no_local_sky'] = bool
        descr['no_local_sky'] = 'If True, turn off local sky model evaluation, but do fit object profile and perform optimal extraction'

        # Masking
        defaults['user_regions'] = None
        dtypes['user_regions'] = [str, list]
        descr['user_regions'] = 'Provides a user-defined mask defining sky regions.  By ' \
                                'default, the sky regions are identified automatically.  To ' \
                                'specify sky regions for *all* slits, provide a comma separated ' \
                                'list of percentages.  For example, setting user_regions = ' \
                                ':10,35:65,80: selects the first 10%, the inner 30%, and the ' \
                                'final 20% of *all* slits as containing sky.  Setting ' \
                                'user_regions = user will attempt to load any SkyRegions ' \
                                'files generated by the user via the pypeit_skysub_regions tool.'

        defaults['mask_by_boxcar'] = False
        dtypes['mask_by_boxcar'] = bool
        descr['mask_by_boxcar'] = 'In global sky evaluation, mask the sky region around the object by the boxcar radius (set in ExtractionPar).'

        defaults['joint_fit'] = False
        dtypes['joint_fit'] = bool
        descr['joint_fit'] = 'Perform a simultaneous joint fit to sky regions using all available slits. ' \
                             'Currently, this parameter is only used for IFU data reduction.'

        defaults['max_mask_frac'] = 0.80
        dtypes['max_mask_frac'] = float
        descr['max_mask_frac'] = 'Maximum fraction of total pixels on a slit that can be masked by the input masks. ' \
                                 'If more than this threshold is masked the code will return zeros and throw a warning.'

        defaults['local_maskwidth'] = 4.0
        dtypes['local_maskwidth'] = float
        descr['local_maskwidth'] = 'Initial width of the region in units of FWHM that will be used for local sky subtraction'


        # Instantiate the parameter set
        super(SkySubPar, self).__init__(list(pars.keys()),
                                        values=list(pars.values()),
                                        defaults=list(defaults.values()),
                                        options=list(options.values()),
                                        dtypes=list(dtypes.values()),
                                        descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])

        # Basic keywords
        parkeys = ['bspline_spacing', 'sky_sigrej', 'global_sky_std', 'no_poly',
                   'user_regions', 'joint_fit', 'mask_by_boxcar',
                   'no_local_sky', 'max_mask_frac', 'local_maskwidth']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for SkySubPar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        pass


class ExtractionPar(ParSet):
    """
    The parameter set used to hold arguments for functionality relevant
    to extraction.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """

    def __init__(self, boxcar_radius=None, std_prof_nsigma=None, sn_gauss=None,
                 model_full_slit=None, skip_extraction=None, skip_optimal=None,
                 use_2dmodel_mask=None, use_user_fwhm=None, return_negative=None):

        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = OrderedDict([(k, values[k]) for k in args[1:]])  # "1:" to skip 'self'

        # Initialize the other used specifications for this parameter
        # set
        defaults = OrderedDict.fromkeys(pars.keys())
        options = OrderedDict.fromkeys(pars.keys())
        dtypes = OrderedDict.fromkeys(pars.keys())
        descr = OrderedDict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set

        # Boxcar Parameters
        defaults['boxcar_radius'] = 1.5
        dtypes['boxcar_radius'] = [int, float]
        descr['boxcar_radius'] = 'Boxcar radius in arcseconds used for boxcar extraction'

        defaults['skip_extraction'] = False
        dtypes['skip_extraction'] = bool
        descr['skip_extraction'] = 'Do not perform an object extraction'

        defaults['skip_optimal'] = False
        dtypes['skip_optimal'] = bool
        descr['skip_optimal'] = 'Perform boxcar extraction only (i.e. skip Optimal and local skysub)'

        defaults['std_prof_nsigma'] = 30.
        dtypes['std_prof_nsigma'] = float
        descr['std_prof_nsigma'] = 'prof_nsigma parameter for Standard star extraction.  Prevents undesired rejection. ' \
                                   'NOTE: Not consumed by the code at present.'

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

        defaults['use_2dmodel_mask'] = True
        dtypes['use_2dmodel_mask'] = bool
        descr['use_2dmodel_mask'] = 'Mask pixels rejected during profile fitting when extracting.' \
                             'Turning this off may help with bright emission lines.'

        defaults['use_user_fwhm'] = False
        dtypes['use_user_fwhm'] = bool
        descr['use_user_fwhm'] = 'Boolean indicating if PypeIt should use the FWHM provided by the user ' \
                                 '(``find_fwhm`` in `FindObjPar`) for the optimal extraction. ' \
                                 'If this parameter is ``False`` (default), PypeIt estimates the FWHM for each ' \
                                 'detected object, and uses ``find_fwhm`` as initial guess.'
        defaults['return_negative'] = False
        dtypes['return_negative'] = bool
        descr['return_negative'] = 'If ``True`` the negative traces will be extracted and saved to disk'

        # Instantiate the parameter set
        super(ExtractionPar, self).__init__(list(pars.keys()),
                                        values=list(pars.values()),
                                        defaults=list(defaults.values()),
                                        options=list(options.values()),
                                        dtypes=list(dtypes.values()),
                                        descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])

        # Basic keywords
        parkeys = ['boxcar_radius', 'std_prof_nsigma', 'sn_gauss', 'model_full_slit',
                   'skip_extraction', 'skip_optimal', 'use_2dmodel_mask', 'use_user_fwhm', 'return_negative']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for ExtractionPar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None

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
    see :ref:`parameters`.
    """
    def __init__(self, calib_dir=None, bpm_usebias=None, biasframe=None, darkframe=None,
                 arcframe=None, tiltframe=None, pixelflatframe=None, pinholeframe=None,
                 alignframe=None, alignment=None, traceframe=None, illumflatframe=None,
                 lampoffflatsframe=None, skyframe=None, standardframe=None, flatfield=None,
                 wavelengths=None, slitedges=None, tilts=None, raise_chk_error=None):


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
        defaults['calib_dir'] = 'Calibrations'
        dtypes['calib_dir'] = str
        descr['calib_dir'] = 'The name of the directory for the processed calibration frames.  ' \
                             'The host path for the directory is set by the redux_path (see ' \
                             ':class:`ReduxPar`).  Beware that success when changing the ' \
                             'default value is not well tested!'

        defaults['raise_chk_error'] = True
        dtypes['raise_chk_error'] = bool
        descr['raise_chk_error'] = 'Raise an error if the calibration check fails'

        defaults['bpm_usebias'] = False
        dtypes['bpm_usebias'] = bool
        descr['bpm_usebias'] = 'Make a bad pixel mask from bias frames? Bias frames must be provided.'

        # Calibration Frames
        defaults['biasframe'] = FrameGroupPar(frametype='bias',
                                              process=ProcessImagesPar(use_biasimage=False,
                                                                       shot_noise=False,
                                                                       use_pixelflat=False,
                                                                       use_illumflat=False,
                                                                       use_specillum=False,
                                                                       combine='median'))
        dtypes['biasframe'] = [ ParSet, dict ]
        descr['biasframe'] = 'The frames and combination rules for the bias correction'

        defaults['darkframe'] = FrameGroupPar(frametype='dark',
                                              process=ProcessImagesPar(use_pixelflat=False,
                                                                       use_illumflat=False,
                                                                       use_specillum=False,
                                                                       mask_cr=True))
        dtypes['darkframe'] = [ ParSet, dict ]
        descr['darkframe'] = 'The frames and combination rules for the dark-current correction'

        # JFH Turning off masking of saturated pixels which causes headaches becauase it was being done unintelligently
        defaults['pixelflatframe'] = FrameGroupPar(frametype='pixelflat',
                                                   process=ProcessImagesPar(satpix='nothing',
                                                                            use_pixelflat=False,
                                                                            use_illumflat=False,
                                                                            use_specillum=False))
        dtypes['pixelflatframe'] = [ ParSet, dict ]
        descr['pixelflatframe'] = 'The frames and combination rules for the pixel flat'

        defaults['illumflatframe'] = FrameGroupPar(frametype='illumflat',
                                                   process=ProcessImagesPar(satpix='nothing',
                                                                            use_pixelflat=False,
                                                                            use_illumflat=False,
                                                                            use_specillum=False))
        dtypes['illumflatframe'] = [ ParSet, dict ]
        descr['illumflatframe'] = 'The frames and combination rules for the illumination flat'

        defaults['lampoffflatsframe'] = FrameGroupPar(frametype='lampoffflats',
                                                     process=ProcessImagesPar(satpix='nothing',
                                                                              use_pixelflat=False,
                                                                              use_illumflat=False,
                                                                              use_specillum=False))
        dtypes['lampoffflatsframe'] = [ ParSet, dict ]
        descr['lampoffflatsframe'] = 'The frames and combination rules for the lamp off flats'

        defaults['pinholeframe'] = FrameGroupPar(frametype='pinhole')
        dtypes['pinholeframe'] = [ ParSet, dict ]
        descr['pinholeframe'] = 'The frames and combination rules for the pinholes'

        defaults['alignframe'] = FrameGroupPar(frametype='align',
                                               process=ProcessImagesPar(satpix='nothing',
                                                                        use_pixelflat=False,
                                                                        use_illumflat=False,
                                                                        use_specillum=False))
        dtypes['alignframe'] = [ ParSet, dict ]
        descr['alignframe'] = 'The frames and combination rules for the align frames'

        defaults['arcframe'] = FrameGroupPar(frametype='arc',
                                             process=ProcessImagesPar(use_pixelflat=False,
                                                                      use_illumflat=False,
                                                                      use_specillum=False))
        dtypes['arcframe'] = [ ParSet, dict ]
        descr['arcframe'] = 'The frames and combination rules for the wavelength calibration'

        defaults['tiltframe'] = FrameGroupPar(frametype='tilt',
                                              process=ProcessImagesPar(use_pixelflat=False,
                                                                       use_illumflat=False,
                                                                       use_specillum=False))
        dtypes['tiltframe'] = [ ParSet, dict ]
        descr['tiltframe'] = 'The frames and combination rules for the wavelength tilts'

        defaults['traceframe'] = FrameGroupPar(frametype='trace',
                                               # Note that CR masking is found to be too problematic!!
                                               process=ProcessImagesPar(use_pixelflat=False,
                                                                        use_illumflat=False,
                                                                        use_specillum=False))

        dtypes['traceframe'] = [ ParSet, dict ]
        descr['traceframe'] = 'The frames and combination rules for images used for slit tracing'

        defaults['standardframe'] = FrameGroupPar(frametype='standard',
                                                  process=ProcessImagesPar(noise_floor=0.01,
                                                                           mask_cr=True))
        dtypes['standardframe'] = [ ParSet, dict ]
        descr['standardframe'] = 'The frames and combination rules for the spectrophotometric ' \
                                 'standard observations'


        defaults['skyframe'] = FrameGroupPar(frametype='sky',
                                                  process=ProcessImagesPar(noise_floor=0.01,
                                                                           mask_cr=True))
        dtypes['skyframe'] = [ ParSet, dict ]
        descr['skyframe'] = 'The frames and combination rules for the sky background ' \
                                 'observations'

        defaults['alignment'] = AlignPar()
        dtypes['alignment'] = [ ParSet, dict ]
        descr['alignment'] = 'Define the procedure for the alignment of traces'

        defaults['flatfield'] = FlatFieldPar()
        dtypes['flatfield'] = [ ParSet, dict ]
        descr['flatfield'] = 'Parameters used to set the flat-field procedure'

        defaults['wavelengths'] = WavelengthSolutionPar()
        dtypes['wavelengths'] = [ ParSet, dict ]
        descr['wavelengths'] = 'Parameters used to derive the wavelength solution'

        defaults['slitedges'] = EdgeTracePar()
        dtypes['slitedges'] = [ ParSet, dict ]
        descr['slitedges'] = 'Slit-edge tracing parameters'

        defaults['tilts'] = WaveTiltsPar()
        dtypes['tilts'] = [ ParSet, dict ]
        descr['tilts'] = 'Define how to trace the slit tilts using the trace frames'

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
        k = np.array([*cfg.keys()])

        # Basic keywords
        parkeys = [ 'calib_dir', 'bpm_usebias', 'raise_chk_error']

        allkeys = parkeys + ['biasframe', 'darkframe', 'arcframe', 'tiltframe', 'pixelflatframe',
                             'illumflatframe', 'lampoffflatsframe',
                             'pinholeframe', 'alignframe', 'alignment', 'traceframe', 'standardframe', 'skyframe',
                             'flatfield', 'wavelengths', 'slitedges', 'tilts']
        badkeys = np.array([pk not in allkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for CalibrationsPar.'.format(k[badkeys]))

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
        pk = 'tiltframe'
        kwargs[pk] = FrameGroupPar.from_dict('tilt', cfg[pk]) if pk in k else None
        pk = 'pixelflatframe'
        kwargs[pk] = FrameGroupPar.from_dict('pixelflat', cfg[pk]) if pk in k else None
        pk = 'illumflatframe'
        kwargs[pk] = FrameGroupPar.from_dict('illumflat', cfg[pk]) if pk in k else None
        pk = 'lampoffflatsframe'
        kwargs[pk] = FrameGroupPar.from_dict('lampoffflats', cfg[pk]) if pk in k else None
        pk = 'pinholeframe'
        kwargs[pk] = FrameGroupPar.from_dict('pinhole', cfg[pk]) if pk in k else None
        pk = 'alignframe'
        kwargs[pk] = FrameGroupPar.from_dict('align', cfg[pk]) if pk in k else None
        pk = 'alignment'
        kwargs[pk] = AlignPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'traceframe'
        kwargs[pk] = FrameGroupPar.from_dict('trace', cfg[pk]) if pk in k else None
        pk = 'standardframe'
        kwargs[pk] = FrameGroupPar.from_dict('standard', cfg[pk]) if pk in k else None
        pk = 'skyframe'
        kwargs[pk] = FrameGroupPar.from_dict('sky', cfg[pk]) if pk in k else None
        pk = 'flatfield'
        kwargs[pk] = FlatFieldPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'wavelengths'
        kwargs[pk] = WavelengthSolutionPar.from_dict(cfg[pk]) if pk in k else None
        pk = 'slitedges'
        kwargs[pk] = EdgeTracePar.from_dict(cfg[pk]) if pk in k else None
        pk = 'tilts'
        kwargs[pk] = WaveTiltsPar.from_dict(cfg[pk]) if pk in k else None

        return cls(**kwargs)

    # TODO: Perform extensive checking that the parameters are valid for
    # the Calibrations class.  May not be necessary because validate will
    # be called for all the sub parameter sets, but this can do higher
    # level checks, if necessary.

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
    see :ref:`parameters`.
    """
    def __init__(self, rdx=None, calibrations=None, scienceframe=None, reduce=None,
                 flexure=None, fluxcalib=None, coadd1d=None, coadd2d=None, sensfunc=None, telluric=None,
                 collate1d=None):

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
        defaults['rdx'] = ReduxPar()
        dtypes['rdx'] = [ ParSet, dict ]
        descr['rdx'] = 'PypeIt reduction rules.'

        defaults['calibrations'] = CalibrationsPar()
        dtypes['calibrations'] = [ ParSet, dict ]
        descr['calibrations'] = 'Parameters for the calibration algorithms'

        defaults['scienceframe'] = FrameGroupPar(frametype='science',
                                                 process=ProcessImagesPar(noise_floor=0.01,
                                                                          mask_cr=True))
        dtypes['scienceframe'] = [ ParSet, dict ]
        descr['scienceframe'] = 'The frames and combination rules for the science observations'

        defaults['reduce'] = ReducePar()
        dtypes['reduce'] = [ParSet, dict]
        descr['reduce'] = 'Parameters determining sky-subtraction, object finding, and ' \
                                'extraction'

        # Flexure is turned OFF by default
        defaults['flexure'] = FlexurePar()
        dtypes['flexure'] = [ParSet, dict]
        descr['flexure'] = 'Parameters used by the flexure-correction procedure.  Flexure ' \
                           'corrections are not performed by default.  To turn on, either ' \
                           'set the parameters in the \'flexure\' parameter group or set ' \
                           '\'flexure = True\' in the \'rdx\' parameter group to use the ' \
                           'default flexure-correction parameters.'

        # Flux calibration is turned OFF by default
        defaults['fluxcalib'] = FluxCalibratePar()
        dtypes['fluxcalib'] = [ParSet, dict]
        descr['fluxcalib'] = 'Parameters used by the flux-calibration procedure.  Flux ' \
                             'calibration is not performed by default.  To turn on, either ' \
                             'set the parameters in the \'fluxcalib\' parameter group or set ' \
                             '\'fluxcalib = True\' in the \'rdx\' parameter group to use the ' \
                             'default flux-calibration parameters.'

        # Coadd1D
        defaults['coadd1d'] = Coadd1DPar()
        dtypes['coadd1d'] = [ParSet, dict]
        descr['coadd1d'] = 'Par set to control 1D coadds.  Only used in the after-burner script.'

        # Coadd2D
        defaults['coadd2d'] = Coadd2DPar()
        dtypes['coadd2d'] = [ParSet, dict]
        descr['coadd2d'] = 'Par set to control 2D coadds.  Only used in the after-burner script.'

        # Sensfunc
        defaults['sensfunc'] = SensFuncPar()
        dtypes['sensfunc'] = [ParSet, dict]
        descr['sensfunc'] = 'Par set to control sensitivity function computation.  Only used in the after-burner script.'

        # Telluric Fit
        defaults['telluric'] = TelluricPar()
        dtypes['telluric'] = [ParSet, dict]
        descr['telluric'] = 'Par set to control telluric fitting.  Only used in the pypeit_sensfunc and pypeit_telluric after-burner scripts.'

        # Collate 1D Fit
        defaults['collate1d'] = Collate1DPar()
        dtypes['collate1d'] = [ParSet, dict]
        descr['collate1d'] = 'Par set to control collating 1d spectra.  Only used in the after-burner script.'

        # Instantiate the parameter set
        super(PypeItPar, self).__init__(list(pars.keys()),
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
                the latter, use `merge_with` to provide one or more
                config files to merge with the defaults to construct the
                full parameter set.
            merge_with (:obj:`str`, :obj:`list`, optional):
                One or more config files with the modifications to
                either default parameters (`cfg_file` is None) or
                the parameters provided by `cfg_file`.  The
                modifications are performed in series so the list order
                of the config files is important.
            evaluate (:obj:`bool`, optional):
                Evaluate the values in the config object before
                assigning them in the subsequent parameter sets.  The
                parameters in the config file are *always* read as
                strings, so this should almost always be true; however,
                see the warning below.

        .. warning::

            When `evaluate` is true, the function runs `eval()` on
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
                (None) and use `merge_with` to provide a set of
                lines to merge with the defaults to construct the full
                parameter set.
            merge_with (:obj:`tuple` or :obj:`list`, optional):
                A tuple containing one more lists
                of strings with lines read, or made to look like
                they are, from a configuration file that should be
                merged with the lines provided by `cfg_lines`, or the
                default parameters.
                The order of the lists in the tuple is important, as
                it sets the order in which the lines are merged.
                Last in line has *highest* priority.
                Or the input may be a list which will be taken 
                as a single item described above.
            evaluate (:obj:`bool`, optional):
                Evaluate the values in the config object before
                assigning them in the subsequent parameter sets.  The
                parameters in the config file are *always* read as
                strings, so this should almost always be true; however,
                see the warning below.

        .. warning::

            When `evaluate` is true, the function runs `eval()` on
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
            # Check it is a tuple
            if isinstance(merge_with, list):
                merge_with = (merge_with,)
            if not isinstance(merge_with, tuple):
                msgs.error('Input merge_with must be a tuple.')
            # Proceed
            for f in merge_with:
                cfg.merge(ConfigObj(f))

        # Evaluate the strings if requested
        if evaluate:
            cfg = util.recursive_dict_evaluate(cfg)

        # Instantiate the object based on the configuration dictionary
        return cls.from_dict(cfg)


    @classmethod
    def from_dict(cls, cfg):
        k = np.array([*cfg.keys()])

        allkeys = ['rdx', 'calibrations', 'scienceframe', 'reduce', 'flexure', 'fluxcalib',
                   'coadd1d', 'coadd2d', 'sensfunc', 'baseprocess', 'telluric', 'collate1d']
        badkeys = np.array([pk not in allkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for PypeItPar.'.format(k[badkeys]))

        kwargs = {}

        pk = 'rdx'
        kwargs[pk] = ReduxPar.from_dict(cfg[pk]) if pk in k else None

        pk = 'calibrations'
        kwargs[pk] = CalibrationsPar.from_dict(cfg[pk]) if pk in k else None

        pk = 'scienceframe'
        kwargs[pk] = FrameGroupPar.from_dict('science', cfg[pk]) if pk in k else None

        pk = 'reduce'
        kwargs[pk] = ReducePar.from_dict(cfg[pk]) if pk in k else None

        # Allow flexure to be turned on using cfg['rdx']
        pk = 'flexure'
        default = FlexurePar()
        kwargs[pk] = FlexurePar.from_dict(cfg[pk]) if pk in k else default

        # Allow flux calibration to be turned on using cfg['rdx']
        pk = 'fluxcalib'
        default = FluxCalibratePar() \
                        if pk in cfg['rdx'].keys() and cfg['rdx']['fluxcalib'] else None
        kwargs[pk] = FluxCalibratePar.from_dict(cfg[pk]) if pk in k else default

        # Allow coadd1d  to be turned on using cfg['rdx']
        pk = 'coadd1d'
        default = Coadd1DPar() \
                        if pk in cfg['rdx'].keys() and cfg['rdx']['coadd1d'] else None
        kwargs[pk] = Coadd1DPar.from_dict(cfg[pk]) if pk in k else default

        # Allow coadd2d  to be turned on using cfg['rdx']
        pk = 'coadd2d'
        default = Coadd2DPar() \
                        if pk in cfg['rdx'].keys() and cfg['rdx']['coadd2d'] else None
        kwargs[pk] = Coadd2DPar.from_dict(cfg[pk]) if pk in k else default

        # Allow coadd2d  to be turned on using cfg['rdx']
        pk = 'sensfunc'
        default = SensFuncPar() \
                        if pk in cfg['rdx'].keys() and cfg['rdx']['sensfunc'] else None
        kwargs[pk] = SensFuncPar.from_dict(cfg[pk]) if pk in k else default

        # Allow telluric to be turned on using cfg['rdx']
        pk = 'telluric'
        default = TelluricPar() \
                        if pk in cfg['rdx'].keys() and cfg['rdx']['telluric'] else None
        kwargs[pk] = TelluricPar.from_dict(cfg[pk]) if pk in k else default

        # collate1d
        pk = 'collate1d'
        default = Collate1DPar()
        kwargs[pk] = Collate1DPar.from_dict(cfg[pk]) if pk in k else default

        if 'baseprocess' not in k:
            return cls(**kwargs)

        # Include any alterations to the basic processing of *all*
        # images
        self = cls(**kwargs)
        baseproc = ProcessImagesPar.from_dict(cfg['baseprocess'])
        self.sync_processing(baseproc)
        return self

    def reset_all_processimages_par(self, **kwargs):
        """
        Change image processing parameter for *all* frame types.

        This function iteratively changes the value of all image processing
        parameters for all frame types in the :class:`CalibrationsPar`, as well
        as the science frames.

        Args:
            **kwargs:
                The list of keywords and values to change for all image
                processing parameters.

        Examples:
            To turn off the slit-illumination correction for all frames:

            >>> from pypeit.spectrographs import load_spectrograph
            >>> spec = load_spectrograph('shane_kast_blue')
            >>> par = spec.default_pypeit_par()
            >>> par.reset_all_processimages_par(use_illumflat=False)

        """
        # Calibrations
        for _key in self['calibrations'].keys():
            if isinstance(self['calibrations'][_key], ParSet) \
                    and 'process' in self['calibrations'][_key].keys():
                for key,value in kwargs.items():
                    self['calibrations'][_key]['process'][key] = value
        # Science frame
        for _key in self.keys():
            if isinstance(self[_key], ParSet) and 'process' in self[_key].keys():
                for key,value in kwargs.items():
                    self[_key]['process'][key] = value

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
    # a full run of PypeIt.  May not be necessary because validate will
    # be called for all the sub parameter sets, but this can do higher
    # level checks, if necessary.
    def validate(self):
        pass

#-----------------------------------------------------------------------------
# Instrument parameters

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
                 diameter=None, eff_aperture=None):

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

        dtypes['eff_aperture'] = [int, float]
        descr['eff_aperture'] = 'Effective aperture of the telescope in m^2'


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
        k = np.array([*cfg.keys()])
        parkeys = [ 'name', 'longitude', 'latitude', 'elevation', 'fratio', 'diameter', 'aperture' ]

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for TelescopePar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    @staticmethod
    def valid_telescopes():
        """
        Return the valid telescopes.
        """
        return [ 'GEMINI-N','GEMINI-S', 'KECK', 'SHANE', 'WHT', 'APF', 'TNG', 'VLT', 'MAGELLAN', 'LBT', 'MMT', 
                'KPNO', 'NOT', 'P200', 'BOK', 'GTC', 'SOAR', 'NTT', 'LDT', 'JWST']

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

    # TODO This method is a place holder until we can get effective apertures for all of the telescopes. I did my best
    # but could not find them all online.
    def eff_aperture(self):
        return np.pi*self['diameter']**2/4.0 if self['eff_aperture'] is None else self['eff_aperture']


class Collate1DPar(ParSet):
    """
    A parameter set holding the arguments for collating, coadding, and archving 1d spectra.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """
    def __init__(self, tolerance=None, dry_run=None, ignore_flux=None, flux=None, match_using=None, exclude_slit_trace_bm=[], exclude_serendip=False, wv_rms_thresh=None, outdir=None, spec1d_outdir=None, refframe=None):

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

        # Threshold for grouping by object
        defaults['tolerance'] = '1.0'
        dtypes['tolerance'] = [str, float]
        descr['tolerance'] = "The tolerance used when comparing the coordinates of objects. If two " \
                             "objects are within this distance from each other, they " \
                             "are considered the same object. If match_using is 'ra/dec' (the default) " \
                             "this is an angular distance. The defaults units are arcseconds but " \
                             "other units supported by astropy.coordinates.Angle can be used " \
                             "(`e.g.`, '0.003d' or '0h1m30s'). If match_using is 'pixel' this is a float."


        # Enables a dry_run in which no coadding is done
        defaults['dry_run'] = False
        dtypes['dry_run'] = bool
        descr['dry_run'] = "If set, the script will display the matching File and Object Ids " \
                           "but will not flux, coadd or archive."

        # Forces coadding using counts instead of flux
        defaults['ignore_flux'] = False
        dtypes['ignore_flux'] = bool
        descr['ignore_flux'] = "If set, the script will only coadd non-fluxed spectra even if flux data is present. " \
                               "Otherwise fluxed spectra are coadded if all spec1ds have been fluxed calibrated."

        # Enables a flux calibration after coadding
        defaults['flux'] = False
        dtypes['flux'] = bool
        descr['flux'] = "If set, the script will flux calibrate using archived sensfuncs before coadding."

        # Directory for output files
        defaults['outdir'] = os.getcwd()
        dtypes['outdir'] = str
        descr['outdir'] = "The path where all coadded output files and report files will be placed."

        # Directory for modified spec1d files
        defaults['spec1d_outdir'] = None
        dtypes['spec1d_outdir'] = str
        descr['spec1d_outdir'] = "The path where all modified spec1d files are placed. These are only created if flux calibration or refframe correction are asked for."

        # What slit flags to exclude
        defaults['exclude_slit_trace_bm'] = []
        dtypes['exclude_slit_trace_bm'] = [list, str]
        descr['exclude_slit_trace_bm'] = "A list of slit trace bitmask bits that should be excluded."

        # What slit flags to exclude
        defaults['exclude_serendip'] = False
        dtypes['exclude_serendip'] = bool
        descr['exclude_serendip'] = "Whether to exclude SERENDIP objects from collating."

        # Wavelength RMS threshold to exclude spectra by.
        defaults['wv_rms_thresh'] = None
        dtypes['wv_rms_thresh'] = float
        descr['wv_rms_thresh'] = "If set, any objects with a wavelength RMS > this value are skipped, else all wavelength RMS values are accepted."

        # How to match objects
        defaults['match_using'] = 'ra/dec'
        options['match_using'] = [ 'pixel', 'ra/dec']
        dtypes['match_using'] = str
        descr['match_using'] = "Determines how 1D spectra are matched as being the same object. Must be either 'pixel' or 'ra/dec'."

        defaults['refframe'] = None
        options['refframe'] = WavelengthSolutionPar.valid_reference_frames()
        dtypes['refframe'] = str
        descr['refframe'] = 'Perform reference frame correction prior to coadding. ' \
                         'Options are: {0}'.format(', '.join(options['refframe']))

        # Instantiate the parameter set
        super(Collate1DPar, self).__init__(list(pars.keys()),
                                           values=list(pars.values()),
                                           defaults=list(defaults.values()),
                                           dtypes=list(dtypes.values()),
                                           descr=list(descr.values()))
        self.validate()

    @classmethod
    def from_dict(cls, cfg):
        k = [*cfg.keys()]
        parkeys = ['tolerance', 'dry_run', 'ignore_flux', 'flux', 'match_using', 'exclude_slit_trace_bm', 'exclude_serendip', 'outdir', 'spec1d_outdir', 'wv_rms_thresh', 'refframe']

        badkeys = np.array([pk not in parkeys for pk in k])
        if np.any(badkeys):
            raise ValueError('{0} not recognized key(s) for Collate1DPar.'.format(k[badkeys]))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = cfg[pk] if pk in k else None
        return cls(**kwargs)

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        pass

