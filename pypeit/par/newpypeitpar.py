"""
Defines parameter sets used to set the behavior for core pypeit
functionality.
"""

from IPython import embed
import numpy as np

from pypeit.par import newparset
from pypeit import log
from pypeit.core import parse


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


class NewFrameGroupPar(newparset.NewParSet):
    """
    New-style parameter set for grouping frames (replacement for FrameGroupPar).
    """
    parameters = {
        'exprng': newparset.set_parameter_definition(
            dtype=list,
            default=[None, None],
            descr=(
                'Used in identifying frames of this type.  This sets the minimum '
                'and maximum allowed exposure times.  There must be two items in '
                'the list.  Use None to indicate no limit; i.e., to select exposures '
                'with any time greater than 30 sec, use exprng = [30, None].'
            ),
        ),
    }

    def validate(self):
        if len(self.data['exprng']) != 2:
            raise ValueError('exprng must be a list with two items.')


class BiasFramePar(NewFrameGroupPar):
    frametype = 'bias'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                use_biasimage=False,
                shot_noise=False,
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
                combine='median',
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class DarkFramePar(NewFrameGroupPar):
    frametype = 'dark'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
                mask_cr=True,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class ScatteredLightFramePar(NewFrameGroupPar):
    frametype = 'scattlight'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                satpix='nothing',
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class PixelFlatFramePar(NewFrameGroupPar):
    frametype = 'pixelflat'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                satpix='nothing',
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class IllumFlatFramePar(NewFrameGroupPar):
    frametype = 'illumflat'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                satpix='nothing',
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class LampOffFlatsFramePar(NewFrameGroupPar):
    frametype = 'lampoffflats'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                satpix='nothing',
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class SlitlessPixFlatBiaFramePar(NewFrameGroupPar):
    frametype = 'slitless_pixflat'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                satpix='nothing',
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
                combine='median',
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class PinholeFramePar(NewFrameGroupPar):
    frametype = 'pinhole'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(),
            descr='Low level parameters used for basic image processing',
        ),
    }


class AlignFramePar(NewFrameGroupPar):
    frametype = 'align'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                satpix='nothing',
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class ArcFramePar(NewFrameGroupPar):
    frametype = 'arc'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class TiltFramePar(NewFrameGroupPar):
    frametype = 'tilt'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class TraceFramePar(NewFrameGroupPar):
    frametype = 'trace'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class StandardFramePar(NewFrameGroupPar):
    frametype = 'standard'
    default_key = f'{frametype}frame'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                noise_floor=0.01,
                mask_cr=True,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class SkyFramePar(NewFrameGroupPar):
    frametype = 'sky'
    default_key = 'skyframe'

    parameters = NewFrameGroupPar.parameters | {
        'process': newparset.set_parameter_definition(
            dtype=NewProcessImagesPar,
            default=NewProcessImagesPar(
                noise_floor=0.01,
                mask_cr=True,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class NewAlignPar(newparset.NewParSet):
    """
    New-style parameter set for alignment tracing (replacement for AlignPar).

    Mirrors the legacy `AlignPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'align'

    parameters = {
        'locations': newparset.set_parameter_definition(
            dtype=[list, np.ndarray],
            default=[0.0, 1.0],
            descr='Locations of the bars, in a list, specified as a fraction of the slit width',
        ),
        'trace_npoly': newparset.set_parameter_definition(
            dtype=int,
            default=4,
            descr='Order of the polynomial to use when fitting the trace of a single bar',
        ),
        'trim_edge': newparset.set_parameter_definition(
            dtype=list,
            default=[0, 0],
            descr='Trim the slit by this number of pixels left/right before finding alignment bars',
        ),
        'snr_thresh': newparset.set_parameter_definition(
            dtype=[int, float],
            default=1.0,
            descr=(
                'S/N ratio threshold for finding an alignment trace. This should be a low '
                'number to ensure that the algorithm finds all bars. The algorithm will '
                'then only use the N most significant detections, where N is the number '
                'of elements specified in the "locations" keyword argument'
            ),
        ),
    }


class NewCoadd1DPar(newparset.NewParSet):
    """
    New-style parameter set for 1D coaddition (replacement for Coadd1DPar).

    Mirrors the legacy `Coadd1DPar` in :mod:`pypeit.par.pypeitpar`.
    """

    valid_extractions = ['BOX', 'OPT']

    valid_wave_methods = ['iref', 'velocity', 'log10', 'linear', 'concatenate']

    valid_scale_methods = ['auto', 'poly', 'median', 'none', 'hand']

    valid_weight_methods = ['auto', 'constant', 'uniform', 'wave_dependent', 'relative', 'ivar']

    default_key = 'coadd1d'

    parameters = {
        'ex_value': newparset.set_parameter_definition(
            dtype=str,
            default='OPT',
            options=valid_extractions,
            descr="The extraction to coadd, i.e. optimal or boxcar. Must be either 'OPT' or 'BOX'",
        ),
        'flux_value': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='If True (default), the code will coadd the fluxed spectra (i.e. the FLAM) in the spec1d files. If False, it will coadd the counts.',
        ),
        'nmaskedge': newparset.set_parameter_definition(
            dtype=int,
            default=2,
            descr='Number of edge pixels to mask. This should be removed/fixed.',
        ),
        'sn_smooth_npix': newparset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Number of pixels to median filter by when computing S/N used to decide how to scale '
                'and weight spectra. If set to None (default), the code will determine the effective '
                'number of good pixels per spectrum in the stack that is being co-added and use 10% of '
                'this neff.'
            ),
        ),
        'sigrej_exp': newparset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Rejection threshold used for rejecting exposures with S/N more than sigrej_exp*sigma '
                'above the median S/N. If None (the default), no rejection is performed. Currently, '
                'only available for multi-slit observations.'
            ),
        ),
        'wave_method': newparset.set_parameter_definition(
            dtype=str,
            default='linear',
            options=valid_wave_methods,
            descr=(
                'Method used to construct wavelength grid for coadding spectra. The routine that creates '
                'the wavelength is :func:`~pypeit.core.wavecal.wvutils.get_wave_grid`. The options are:'
                " "
                "'iref' -- Use the first wavelength array.  "
                "'velocity' -- Grid is uniform in velocity.  "
                "'log10' -- Grid is uniform in log10(wave). This is the same as velocity.  "
                "'linear' -- Grid is uniform in lambda.  "
                "'concatenate' -- Meld the input wavelength arrays"
            ),
        ),
        'dv': newparset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                "Dispersion in units of km/s in case you want to specify it in the get_wave_grid  (for the 'velocity' option), "
                "otherwise a median value is computed from the data."
            ),
        ),
        'dwave': newparset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                "Dispersion in Angstroms in case you want to specify it in the get_wave_grid  (for the 'linear' option), "
                "otherwise a median value is computed from the data."
            ),
        ),
        'dloglam': newparset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                "Dispersion in units of log10(wave) in case you want to specify it in the get_wave_grid  (for the 'velocity' or 'log10' options), "
                "otherwise a median value is computed from the data."
            ),
        ),
        'wave_grid_min': newparset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Used in case you want to specify the minimum wavelength in your wavelength grid, default=None computes from data'
            ),
        ),
        'wave_grid_max': newparset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Used in case you want to specify the maximum wavelength in your wavelength grid, default=None computes from data'
            ),
        ),
        'spec_samp_fact': newparset.set_parameter_definition(
            dtype=float,
            default=1.0,
            descr=(
                "Make the wavelength grid  sampling finer (spec_samp_fact < 1.0) or coarser "
                "(spec_samp_fact > 1.0) by this sampling factor. This basically multiples the 'native' "
                "spectral pixels by spec_samp_fact, i.e. units spec_samp_fact are pixels."
            ),
        ),
        'ref_percentile': newparset.set_parameter_definition(
            dtype=[int, float],
            default=70.0,
            descr=(
                'Percentile used for selecting the minimum SNR cut from a reference spectrum used to '
                'robustly determine the median ratio between spectra. This parameter is used by '
                'coadd1d.robust_median_ratio as part of the automatic rescaling procedure. Pixels '
                'above this percentile cut are deemed the "good" pixels and are used to compute the '
                'ratio of two spectra.  This must be a number between 0 and 100.'
            ),
        ),
        'maxiter_scale': newparset.set_parameter_definition(
            dtype=int,
            default=5,
            descr='Maximum number of iterations performed for rescaling spectra.',
        ),
        'sigrej_scale': newparset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr='Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec.',
        ),
        'scale_method': newparset.set_parameter_definition(
            dtype=str,
            default='auto',
            options=valid_scale_methods,
            descr=(
                "Method used to rescale the spectra prior to coadding. The options are:" 
                " "
                "'auto' -- Determine the scaling method automatically based on the S/N ratio which works well.  "
                "'poly' -- Polynomial rescaling.  " 
                "'median' -- Median rescaling  " 
                "'none' -- Do not rescale.  " 
                "'hand' -- Pass in hand scaling factors. This option is not well tested."
            ),
        ),
        'sn_min_medscale': newparset.set_parameter_definition(
            dtype=[int, float],
            default=0.5,
            descr='For scale method set to ``auto``, this sets the minimum SNR for which median scaling is attempted.',
        ),
        'sn_min_polyscale': newparset.set_parameter_definition(
            dtype=[int, float],
            default=2.0,
            descr='For scale method set to ``auto``, this sets the minimum SNR for which polynomial scaling is attempted.',
        ),
        'weight_method': newparset.set_parameter_definition(
            dtype=str,
            default='auto',
            options=valid_weight_methods,
            descr=(
                "Method used to weight the spectra for coadding. The options are:" 
                " " 
                "'auto' -- Use constant weights if rms_sn < 3.0, otherwise use wavelength dependent." 
                "'constant' -- Constant weights based on rms_sn**2" 
                "'uniform' --  Uniform weighting" 
                "'wave_dependent' -- Wavelength dependent weights will be used irrespective of the rms_" 
                "sn ratio. This option will not work well at low S/N ratio although it is useful for " 
                "objects where only a small fraction of the spectral coverage has high S/N ratio " 
                "(like high-z quasars)." 
                "'relative' -- Apply relative weights implying one reference exposure will receive unit " 
                "weight at all wavelengths and all others receive relatively wavelength dependent "
                "weights . Note, relative weighting will only work well " 
                "when there is at least one spectrum with a reasonable S/N, and a continuum. " 
                "This option may only be better when the object being used has a strong " 
                "continuum + emission lines. This is particularly useful if you " 
                "are dealing with highly variable spectra (e.g. emission lines) and" 
                "require a precision better than ~1 per cent." 
                "'ivar' -- Use inverse variance weighting. This is not well tested and should probably be deprecated."
            ),
        ),
        'maxiter_reject': newparset.set_parameter_definition(
            dtype=int,
            default=5,
            descr='Maximum number of iterations for stacking and rejection. The code stops iterating ' \
                  'either when the output mask does not change betweeen successive iterations or when ' \
                  'maxiter_reject is reached.',
        ),
        'lower': newparset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr='Lower rejection threshold used for rejecting pixels when combining spectra in units of sigma.',
        ),
        'upper': newparset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr='Upper rejection threshold used for rejecting pixels when combining spectra in units of sigma.',
        ),
        'maxrej': newparset.set_parameter_definition(
            dtype=int,
            descr=(
                'Coadding performs iterative rejection by comparing each exposure to a preliminary stack of '
                'all the exposures. If this parameter is set then it will not reject more than maxrej pixels '
                'per iteration of this rejection. The default is None, which means no maximum on rejected pixels.'
            ),
        ),
        'sn_clip': newparset.set_parameter_definition(
            dtype=[int, float],
            default=30.0,
            descr=(
                'Errors are capped during rejection so that the S/N is never greater than sn_clip. This '
                'prevents overly aggressive rejection in high S/N ratio spectrum which neverthless differ '
                'at a level greater than the formal S/N due to systematics.'
            ),
        ),
        'nbests': newparset.set_parameter_definition(
            dtype=[list, int],
            descr=(
                'Number of orders to use for estimating the per exposure weights. Default is None, '
                'which will just use one fourth of the total number of orders. This is only used for Echelle'
            ),
        ),
        'filter': newparset.set_parameter_definition(
            dtype=str,
            default='none',
            descr='Filter for scaling.  See flux_calib.load_fitler_file() for naming.  Ignore if none',
        ),
        'mag_type': newparset.set_parameter_definition(
            dtype=str,
            default='AB',
            descr='Magnitude type.  AB is the only option currently allowed',
        ),
        'filter_mag': newparset.set_parameter_definition(
            dtype=float,
            descr='Magnitude of the source in the given filter',
        ),
        'filter_mask': newparset.set_parameter_definition(
            dtype=[str, list],
            descr=(
                'List of wavelength regions to mask when doing the scaling (`i.e.`, occasional junk pixels). '
                'Colon and comma separateed, e.g.   5552:5559,6010:6030'
            ),
        ),
        'coaddfile': newparset.set_parameter_definition(
            dtype=str,
            descr='Output filename',
        ),
    }

    def validate(self):
        allowed_extensions = self.valid_ex()
        if self.data['ex_value'] not in allowed_extensions:
            raise ValueError("'ex_value' must be one of:\n" + ", ".join(allowed_extensions))

        allowed_wave_methods = self.valid_wave_methods()
        if self.data['wave_method'] not in allowed_wave_methods:
            raise ValueError("'wave_method' must be one of:\n" + ", ".join(allowed_wave_methods))

        allowed_scale_methods = self.valid_scale_methods()
        if self.data['scale_method'] not in allowed_scale_methods:
            raise ValueError("'scale_method' must be one of:\n" + ", ".join(allowed_scale_methods))

        allowed_weight_methods = self.valid_weight_methods()
        if self.data['weight_method'] not in allowed_weight_methods:
            raise ValueError("'weight_method' must be one of:\n" + ", ".join(allowed_weight_methods))


class NewCoadd2DPar(newparset.NewParSet):
    """
    New-style parameter set for 2D coaddition (replacement for Coadd2DPar).

    Mirrors the legacy `Coadd2DPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'coadd2d'

    parameters = {
        'only_slits': newparset.set_parameter_definition(
            dtype=[str, list],
            default=None,
            descr=(
                'Restrict coaddition to one or more of slits. Example syntax -- '
                'DET01:175,DET02:205 or MSC02:2234. This and ``exclude_slits`` '
                'are mutually exclusive. If both are provided, ``only_slits`` takes precedence.'
            ),
        ),
        'exclude_slits': newparset.set_parameter_definition(
            dtype=[str, list],
            default=None,
            descr=(
                'Exclude one or more slits from the coaddition. Example syntax -- '
                'DET01:175,DET02:205 or MSC02:2234. This and ``only_slits`` '
                'are mutually exclusive. If both are provided, ``only_slits`` takes precedence.'
            ),
        ),
        'offsets': newparset.set_parameter_definition(
            dtype=[str, list],
            default='auto',
            descr=(
                'Offsets for the images being combined (spat pixels). Options are: '
                '``maskdef_offsets``, ``header``, ``auto``, and a list of offsets. '
                'Use ``maskdef_offsets`` to use the offsets computed during the slitmask design matching '
                '(currently available for these :ref:`slitmask_info_instruments` only). If equal '
                'to ``header``, the dither offsets recorded in the header, when available, will be used. '
                'If ``auto`` is chosen, PypeIt will try to compute the offsets using a reference object '
                'with the highest S/N, or using a list of object ids selected by the user (see ``user_obj_ids``). '
                'If a list of offsets is provided, PypeIt will use it.'
            ),
        ),
        'spat_toler': newparset.set_parameter_definition(
            dtype=int,
            default=5,
            descr='This parameter provides the desired tolerance in spatial pixel used ' \
                  'to identify slits in different exposures',
        ),
        'use_slits4wvgrid': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='If True, use the slits to set the trace down the center',
        ),
        'weights': newparset.set_parameter_definition(
            dtype=[str, list],
            default='auto',
            descr=(
                'Mode for the weights used to coadd images. Options are: ' 
                '``auto``, ``uniform``, or a list of weights. ' 
                'If a list of weights is provided, PypeIt will use it.' 
                'if ``uniform`` is used, uniform weights will be applied.' 
                'If ``auto`` is used, PypeIt will try to compute the weights ' 
                'using a reference object with the highest S/N, or using a list ' 
                'of object ids selected by the user (see ``user_obj_ids``). If the reference object ' 
                'is not found, the code will use uniform weights. '
            ),
        ),
        'user_obj_ids': newparset.set_parameter_definition(
            dtype=list,
            default=None,
            descr=(
                'List of unique object identifiers that the user wants to use '
                'to compute the weights and/or the offsets for coadding images. '
                'For longslit/multislit spectroscopy, provide the ``SPAT_PIXPOS_ID`` '
                'of the object in each of the exposures. For echelle spectroscopy, '
                'provide the ``ECH_FRACPOS_ID`` of the object in each exposure. ' \
                'These unique object identifiers can be found in the spec1d*.txt ' \
                'files for each exposure. See :doc:`out_spec1D` for more info about ' \
                '``SPAT_PIXPOS_ID`` and ``ECH_FRACPOS_ID``. This parameter must always ' \
                'be a list of the same length as the number of exposures being coadded. ' \
                'If this parameter is not ``None``, it will be used to compute the offsets ' \
                'only if ``offsets = auto``, and it will used to compute the weights ' \
                'only if ``weights = auto``.'
            ),
        ),
        'manual': newparset.set_parameter_definition(
            dtype=[str, list],
            default=None,
            descr=(
                'Manual extraction parameters for eac aperture to extract.  For a ' 
                'single detector, use det:spat:spec:fwhm:boxcar_radius.  For a ' 
                'mosiac, use ``(det1,det2,...):spat:spec:fwhm:boxcar_radius``, where ' 
                '``(det1,det2,...)`` is the list of detectors in the mosaic.  ' 
                'Multiple manual extraction apertures are separated by semicolons; ' 
                'e.g., ``(1,2,3):22.4:608.1:3.; (1,2,3):82.4:608.1:3.``.  Note ' 
                '``spat,spec`` are in the pixel coordinates of the pseudo-image ' 
                'generated by COADD2D; ``fwhm`` is in pixels, and ``boxcar_radius`` ' 
                'is optional and **in pixels (not arcsec!)**.'
            ),
        ),
        'wave_method': newparset.set_parameter_definition(
            dtype=str,
            default=None,
            options=['iref', 'velocity', 'log10', 'linear'],
            descr=(
                "Argument to :func:`~pypeit.core.wavecal.wvutils.get_wave_grid` method, which determines how "
                "the 2d coadd wavelength grid is constructed. The default is None, which will use a linear grid" 
                "for longslit/multislit coadds and a log10 grid for echelle coadds. " 
                "Currently supported options with 2d coadding are:" 
                "* 'iref' -- Use one of the exposures (the first) as the reference for the wavelength grid " 
                "* 'velocity' -- Grid is uniform in velocity" 
                "* 'log10'  -- Grid is uniform in log10(wave). This is the same as velocity." 
                "* 'linear' -- Grid is uniform in wavelength"
            ),
        ),
        'spec_samp_fact': newparset.set_parameter_definition(
            dtype=float,
            default=1.0,
            descr=(
                "Make the wavelength grid sampling finer (``spec_samp_fact`` less than 1.0)" \
                "or coarser (``spec_samp_fact`` greater than 1.0) by this sampling factor." \
                "This  multiples the 'native' spectral pixel size by ``spec_samp_fact``," \
                "i.e. the units of ``spec_samp_fact`` are pixels."
            ),
        ),
        'spat_samp_fact': newparset.set_parameter_definition(
            dtype=float,
            default=1.0,
            descr=(
                "Make the spatial sampling finer (``spat_samp_fact`` less" \
                "than 1.0) or coarser (``spat_samp_fact`` greather than 1.0) by" \
                "this sampling factor. This basically multiples the 'native'" \
                "spatial pixel size by ``spat_samp_fact``, i.e. the units of" \
                "``spat_samp_fact`` are pixels."
            ),
        ),
    }

    def validate(self):
        allowed_wave_methods = ['iref', 'velocity', 'log10', 'linear']
        if self.data['wave_method'] is not None and self.data['wave_method'] not in allowed_wave_methods:
            raise ValueError(f"If 'wave_method' is set it must be one of: {', '.join(allowed_wave_methods)}")

        # Normalize manual extraction entries if provided
        if self.data['manual'] is not None:
            self.data['manual'] = ';'.join(parse.fix_config_par_image_location(self.data['manual']))

