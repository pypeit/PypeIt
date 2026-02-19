"""
Defines parameter sets used to set the behavior for core pypeit
functionality.
"""

from IPython import embed

from pypeit.par import newparset
from pypeit import log


class NewScatteredLightPar(newparset.NewParSet):
    """
    The parameter set used to hold arguments for modelling the scattered light.

    For a table with the current keywords, defaults, and descriptions,
    see :ref:`parameters`.
    """

    valid_scattlight_methods = ['model', 'frame', 'archive']
    """
    Return the valid scattered light methods.
    """

    valid_finecorr_scattlight_methods = ['median', 'poly']
    """
    Return the valid scattered light methods.
    """

    default_key = 'scattlight'

    card_prefix = 'SCLT'

    parameters = {
        'method': newparset.set_parameter_definition(
            dtype=str,
            default='model',
            options=valid_scattlight_methods,
            descr=(
                'Method used to fit the overscan.  Options are '
                f'{", ".join(valid_scattlight_methods)}.  Select "model" '
                'to use the scattered light model parameters derived from a user-specified frame '
                'during the reduction (you will need to identify appropriate "scattlight" frames '
                'in your pypeit file).  Select "frame" to use each individual frame to determine '
                'the scattered light directly from that frame.  Select "archive" to use an '
                'archival model parameter solution for the scattered light (this option is not '
                'currently available for all spectrographs).'
            )
        ),
        'finecorr_method': newparset.set_parameter_definition(
            dtype=str,
            options=valid_finecorr_scattlight_methods,
            descr=(
                'Method used for the "fine correction" for scattered light.  This can be None to '
                'skip the correction; the other options are '
                f'{", ".join(valid_finecorr_scattlight_methods)}.  Select "median" to subtract a '
                'constant value from an entire CCD row, based on a median of the pixels that are '
                'not on slits (see also, "finecorr_pad").  Select "poly" to fit a polynomial to '
                'the scattered light in each row, based on the pixels that are not on slits '
                '(see also, "finecorr_pad").'
            )
        ),
        'finecorr_pad': newparset.set_parameter_definition(
            dtype=int,
            default=4,
            descr=(
                'Number of unbinned pixels by which to extend the slit edges by when masking the '
                'slits for the fine correction to the scattered light.'
            )
        ),
        'finecorr_order': newparset.set_parameter_definition(
            dtype=int,
            default=2,
            descr=(
                'Polynomial order to use for the fine correction to the scattered light '
                'subtraction. It should be a low value.'
            )
        ),
        'finecorr_mask': newparset.set_parameter_definition(
            dtype=[int, list],
            descr=(
                'The inter-slit regions to mask during the fine correction to the scattered '
                'light.  Each integer corresponds to an inter-slit region.  For example, "0" '
                'corresponds to all pixels left of the leftmost slit, whereas "1" corresponds to '
                'all pixels between the first and second slit (counting from the left).  Provide '
                'either a single integer value or a list of integer values. The default (None) '
                'means that no inter-slit regions will be masked.'
            )
        ),
    }


class NewProcessImagesPar(newparset.NewParSet):
    """
    New-style parameter set for basic image processing using `newparset.NewParSet`.

    This replaces the old instance-driven __init__ with a class-level
    `parameters` specification. The `newparset.NewParSet` base class handles defaulting,
    type/options validation, and instantiation.
    """

    valid_overscan_methods = ['chebyshev', 'polynomial', 'savgol', 'median', 'odd_even']
    """
    Valid methods for the overscan subtraction.
    """

    valid_combine_methods = ['median', 'mean']
    """
    Valid methods for image combining.
    """

    valid_saturation_handling = ['reject', 'force', 'nothing']
    """
    Valid methods for dealing with saturated pixels.
    """

    default_key = 'process'

    card_prefix = 'IPRC'

    parameters = {
        'trim': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Trim the image to the detector supplied region',
        ),
        'apply_gain': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Convert the ADUs to electrons using the detector gain',
        ),
        'orient': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Orient the raw image into the PypeIt frame',
        ),
        'use_biasimage': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Use a bias image.  If True, one or more must be supplied in the PypeIt file.',
        ),
        'use_overscan': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Subtract off the overscan.  Detector *must* have one or code will crash.',
        ),
        'overscan_method': newparset.set_parameter_definition(
            dtype=str,
            default='savgol',
            options=valid_overscan_methods,
            descr=(
                'Method used to fit the overscan. '
                f'Options are: {", ".join(valid_overscan_methods)}  Note: Method "polynomial" '
                'is identical to "chebyshev"; the former is deprecated and will be removed.'
            ),
        ),
        'overscan_par': newparset.set_parameter_definition(
            dtype=[int, list],
            default=[5, 65],
            descr=(
                'Parameters for the overscan subtraction.  For '
                "'chebyshev' or 'polynomial', set overcan_par = order; "
                "for 'savgol', set overscan_par = order, window size ; "
                "for 'median', set overscan_par = None or omit the keyword."
            ),
        ),
        'correct_nonlinear': newparset.set_parameter_definition(
            dtype=list,
            default=None,
            descr=(
                'Correct for non-linear response of the detector.  If None, '
                'no correction is performed. If a list, then the list should be '
                'the non-linear correction parameter (alpha), where the functional '
                'form is given by Ct = Cm (1 + alpha x Cm), with Ct and Cm the true '
                'and measured counts. This parameter is usually '
                'hard-coded for a given spectrograph, and should otherwise be left as None.'
            ),
        ),
        'use_darkimage': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Subtract off a dark image.  If True, one or more darks must be provided.',
        ),
        'dark_expscale': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'If designated dark frames are used and have a different '
                'exposure time than the science frames, scale the counts by '
                'the by the ratio in the exposure times to adjust the dark '
                'counts for the difference in exposure time.  WARNING: You '
                'should always take dark frames that have the same exposure '
                'time as your science frames, so use this option with care!'
            ),
        ),
        'use_pattern': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Subtract off a detector pattern. This pattern is assumed to be '
                'sinusoidal along one direction, with a frequency that is '
                'constant across the detector.'
            ),
        ),
        'subtract_continuum': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Subtract off the continuum level from an image. This parameter should only '
                'be set to True to combine arcs with multiple different lamps. '
                'For all other cases, this parameter should probably be False.'
            ),
        ),
        'subtract_scattlight': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Subtract off the scattered light from an image. This parameter should only '
                'be set to True for spectrographs that have dedicated methods to subtract '
                'scattered light. For all other cases, this parameter should be False.'
            ),
        ),
        'scattlight': newparset.set_parameter_definition(
            dtype=NewScatteredLightPar,
            default=NewScatteredLightPar(),
            descr='Scattered light subtraction parameters.',
        ),
        'empirical_rn': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'If True, use the standard deviation in the overscan region to '
                'measure an empirical readnoise to use in the noise model.'
            ),
        ),
        'shot_noise': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Use the bias- and dark-subtracted image to calculate and include '
                'electron count shot noise in the image processing error budget'
            ),
        ),
        'noise_floor': newparset.set_parameter_definition(
            dtype=float,
            default=0.0,
            descr=(
                'Impose a noise floor by adding the provided fraction of the '
                'bias- and dark-subtracted electron counts to the error budget.  '
                'E.g., a value of 0.01 means that the S/N of the counts in the '
                'image will never be greater than 100.'
            ),
        ),
        'use_pixelflat': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Use the pixel flat to make pixel-level corrections.  A '
                'pixelflat image must be provied.'
            ),
        ),
        'use_illumflat': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Use the illumination flat to correct for the illumination profile of each slit.',
        ),
        'use_specillum': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Use the relative spectral illumination profiles to correct '
                'the spectral illumination profile of each slit. This is '
                'primarily used for slicer IFUs.  To use this, you must set '
                '``slit_illum_relative=True`` in the ``flatfield`` parameter set!'
            ),
        ),
        'skip_write_2d': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Skip writing the 2D spectrum for science frames.  WARNING: '
                'This option should only be considered for reducing the volume '
                'of output data when processing large numbers of frames and only '
                'after ensuring the quality of the resulting reductions.'
            ),
        ),
        'spat_flexure_correct': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Correct slits, illumination flat, etc. for flexure',
        ),
        'spat_flexure_maxlag': newparset.set_parameter_definition(
            dtype=int,
            default=20,
            descr='Maximum of possible spatial flexure correction, in pixels',
        ),
        'spat_flexure_sigdetect': newparset.set_parameter_definition(
            dtype=[int, float],
            default=5.0,
            descr=(
                'Sigma threshold above fluctuations in the '
                'Sobel-filtered significance image, used for '
                'finding slit edges in the spectral image, '
                'for which the spatial flexure is computed.'
            ),
        ),
        'spat_flexure_vrange': newparset.set_parameter_definition(
            dtype=tuple,
            default=None,
            descr=(
                'This parameter is used when generating the QA plot for the spatial flexure. '
                'It sets the data range (vmin,vmax) used by the colormap when showing the '
                'spectral image. If None, the range is set automatically.'
            ),
        ),
        'combine': newparset.set_parameter_definition(
            dtype=str,
            default='mean',
            options=valid_combine_methods,
            descr=(
                'Method used to combine multiple frames.  Options are: '
                f'{", ".join(valid_combine_methods)}'
            ),
        ),
        'clip': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Perform sigma clipping when combining.  Only used with combine=mean',
        ),
        'scale_to_mean': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='If True, scale the input images to have the same mean before combining.',
        ),
        'comb_sigrej': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Sigma-clipping level for when clip=True; '
                'Use None for automatic limit (recommended).  '
            ),
        ),
        'satpix': newparset.set_parameter_definition(
            dtype=str,
            default='reject',
            options=valid_saturation_handling,
            descr=(
                'Handling of saturated pixels.  Options are: '
                f'{", ".join(valid_saturation_handling)}'
            ),
        ),
        'mask_cr': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Identify CRs and mask them',
        ),
        'n_lohi': newparset.set_parameter_definition(
            dtype=list,
            default=[0, 0],
            descr=(
                'Number of pixels to reject at the lowest and highest ends of the '
                'distribution; i.e., n_lohi = low, high.  Use None for no limit.'
            ),
        ),
        'lamaxiter': newparset.set_parameter_definition(
            dtype=int,
            default=1,
            descr='Maximum number of iterations for LA cosmics routine.',
        ),
        'grow': newparset.set_parameter_definition(
            dtype=[int, float],
            default=1.5,
            descr=(
                'Factor by which to expand regions with cosmic rays detected by the '
                'LA cosmics routine.'
            ),
        ),
        'rmcompact': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Remove compact detections in LA cosmics routine',
        ),
        'sigclip': newparset.set_parameter_definition(
            dtype=[int, float],
            default=4.5,
            descr='Sigma level for rejection in LA cosmics routine',
        ),
        'sigfrac': newparset.set_parameter_definition(
            dtype=[int, float],
            default=0.3,
            descr='Fraction for the lower clipping threshold in LA cosmics routine.',
        ),
        'objlim': newparset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr='Object detection limit in LA cosmics routine',
        ),
    }

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['n_lohi'] is not None and len(self.data['n_lohi']) != 2:
            raise ValueError('n_lohi must be a list of two numbers.')

        # Check the overscan methods
        if self.data['use_overscan']:
            if self.data['overscan_par'] is None:
                raise ValueError('No overscan method parameters defined!')

            # Convert param to list
            # TODO: We should impose this using the dtype
            if isinstance(self.data['overscan_par'], int):
                self.data['overscan_par'] = [self.data['overscan_par']]

            if (
                self.data['overscan_method'] in ['polynomial', 'chebyshev']
                and len(self.data['overscan_par']) != 1
            ):
                raise ValueError(
                    'For chebyshev/polynomial overscan method, set overscan_par = order'
                )

            if self.data['overscan_method'] == 'savgol' and len(self.data['overscan_par']) != 2:
                raise ValueError(
                    'For savgol overscan method, set overscan_par = order, window size'
                )

            if self.data['overscan_method'] == 'median' and self.data['overscan_par'] is not None:
                log.warning('No parameters necessary for median overscan method.  Ignoring input.')

        # Check the consistency of the flat-fielding approach
        if (
            not self.data['use_pixelflat']
            and (self.data['use_illumflat'] or self.data['use_specillum'])
        ):
            raise ValueError(
                'To apply a slit-illumination or spectral flat-field correction, you must also '
                'apply the pixel-flat correction.'
            )

