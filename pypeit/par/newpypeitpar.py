"""
Defines parameter sets used to set the behavior for core pypeit
functionality.
"""
from pathlib import Path

from IPython import embed
import numpy as np
import os

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


class SlitlessPixFlatFramePar(NewFrameGroupPar):
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


class ScienceFramePar(NewFrameGroupPar):
    frametype = 'science'
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
                'of object ids selected by the user indicating a reference object '
                'in each exposure (see ``user_obj_ids``). If the reference object '
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
        # Normalize manual extraction entries if provided
        if self.data['manual'] is not None:
            self.data['manual'] = ';'.join(parse.fix_config_par_image_location(self.data['manual']))


class NewCubePar(newparset.NewParSet):
    """
    New-style parameter set for cube generation (replacement for CubePar).

    Mirrors the legacy `CubePar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'cube'

    valid_weight_methods = ['auto', 'constant', 'uniform', 'wave_dependent', 'relative', 'ivar']

    parameters = {
        'slit_spec': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If the data use slits in one spatial direction, set this to True. '
                'If the data uses fibres for all spaxels, set this to False.'
            ),
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
        'align': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'If set to True, the input frames will be spatially aligned by cross-correlating the '
                'whitelight images with either a reference image (see ``reference_image``) or the whitelight '
                'image that is generated using the first spec2d listed in the coadd3d file. Alternatively, '
                'the user can specify the offsets (i.e. Delta RA x cos(dec) and Delta Dec, both in arcsec) '
                'in the spec2d block of the coadd3d file. See the documentation for examples of this usage.'
            ),
        ),
        'combine': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='If set to True, the input frames will be combined. Otherwise, a separate datacube will be generated for each input spec2d file, and will be saved as a spec3d file.',
        ),
        'output_filename': newparset.set_parameter_definition(
            dtype=str,
            default="",
            descr=(
                'If combining multiple frames, this string sets the output filename of '
                'the combined datacube. If combine=False, the output filenames will be '
                'prefixed with ``spec3d_*``'
            ),
        ),
        'sensfile': newparset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=(
                'Filename of a sensitivity function to use to flux calibrate your datacube. '
                'The sensitivity function file will also be used to correct the relative scales '
                'of the slits.'
            ),
        ),
        'reference_image': newparset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=(
                'White light image of a previously combined datacube. The white light '
                'image will be used as a reference when calculating the offsets of the '
                'input spec2d files. Ideally, the reference image should have the same '
                'shape as the data to be combined (i.e. set the ra_min, ra_max etc. params '
                'so they are identical to the reference image).'
            ),
        ),
        'save_whitelight': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Save a white light image of the combined datacube. The output filename '
                'will be given by the "output_filename" variable with a suffix "_whitelight". '
                'Note that the white light image collapses the flux along the wavelength axis, '
                'so some spaxels in the 2D white light image may have different wavelength '
                'ranges. To set the wavelength range, use the "whitelight_range" parameter. '
                'If combine=False, the individual spec3d files will have a suffix "_whitelight".'
            ),
        ),
        'whitelight_range': newparset.set_parameter_definition(
            dtype=list,
            default=[None, None],
            descr=(
                'A two element list specifying the wavelength range over which to generate the '
                'white light image. The first (second) element is the minimum (maximum) '
                'wavelength to use. If either of these elements are None, PypeIt will '
                'automatically use a wavelength range that ensures all spaxels have the '
                'same wavelength coverage. Note, if you are using a reference_image to align '
                'all frames, it is preferable to use the same white light wavelength range '
                'for all white light images. For example, you may wish to use an emission '
                'line map to register two frames.'
            ),
        ),
        'method': newparset.set_parameter_definition(
            dtype=str,
            default='subpixel',
            options=['subpixel', 'ngp'],
            descr=(
                'What method should be used to generate the datacube. There are currently two options: '
                '(1) "subpixel" (default) - this algorithm divides each pixel in the spec2d frames '
                'into subpixels, and assigns each subpixel to a voxel of the datacube. Flux is conserved, '
                'but voxels are correlated, and the error spectrum does not account for covariance between '
                'adjacent voxels. See also, spec_subpixel and spat_subpixel. '
                '(2) "ngp" (nearest grid point) - this algorithm is effectively a 3D histogram. Flux is '
                'conserved, voxels are not correlated, however this option suffers the same downsides as '
                'any histogram; the choice of bin sizes can change how the datacube appears. This algorithm '
                'takes each pixel on the spec2d frame and puts the flux of this pixel into one voxel in the '
                'datacube. Depending on the binning used, some voxels may be empty (zero flux) while a '
                'neighboring voxel might contain the flux from two spec2d pixels. Note that all spec2d '
                'pixels that contribute to the same voxel are inverse variance weighted (e.g. if two '
                'pixels have the same variance, the voxel would be assigned the average flux of the two '
                'pixels).'
            ),
        ),
        'spec_subpixel': newparset.set_parameter_definition(
            dtype=int,
            default=5,
            descr=(
                'When method=subpixel, spec_subpixel sets the subpixellation scale of '
                'each detector pixel in the spectral direction. The total number of subpixels '
                'in each pixel is given by spec_subpixel x spat_subpixel. The default option '
                'is to divide each spec2d pixel into 25 subpixels during datacube creation. '
                'See also, spat_subpixel and slice_subpixel.'
            ),
        ),
        'spat_subpixel': newparset.set_parameter_definition(
            dtype=int,
            default=5,
            descr=(
                'When method=subpixel, spat_subpixel sets the subpixellation scale of '
                'each detector pixel in the spatial direction. The total number of subpixels '
                'in each pixel is given by spec_subpixel x spat_subpixel. The default option '
                'is to divide each spec2d pixel into 25 subpixels during datacube creation. '
                'See also, spec_subpixel and slice_subpixel.'
            ),
        ),
        'slice_subpixel': newparset.set_parameter_definition(
            dtype=int,
            default=5,
            descr=(
                'When method=subpixel, slice_subpixel sets the subpixellation scale of '
                'each IFU slice. The default option is to divide each slice into 5 sub-slices '
                'during datacube creation. See also, spec_subpixel and spat_subpixel.'
            ),
        ),
        'ra_min': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Minimum RA to use when generating the WCS. If None, the default is minimum RA '
                'based on the WCS of all spaxels. Units should be degrees.'
            ),
        ),
        'ra_max': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Maximum RA to use when generating the WCS. If None, the default is maximum RA '
                'based on the WCS of all spaxels. Units should be degrees.'
            ),
        ),
        'dec_min': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Minimum DEC to use when generating the WCS. If None, the default is minimum DEC '
                'based on the WCS of all spaxels. Units should be degrees.'
            ),
        ),
        'dec_max': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Maximum DEC to use when generating the WCS. If None, the default is maximum DEC '
                'based on the WCS of all spaxels. Units should be degrees.'
            ),
        ),
        'wave_min': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Minimum wavelength to use when generating the WCS. If None, the default is '
                'minimum wavelength based on the WCS of all spaxels. Units should be Angstroms.'
            ),
        ),
        'wave_max': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Maximum wavelength to use when generating the WCS. If None, the default is '
                'maximum wavelength based on the WCS of all spaxels. Units should be Angstroms.'
            ),
        ),
        'spatial_delta': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'The spatial size of each spaxel to use when generating the WCS (in arcsec). '
                'If None, the default is set by the spectrograph file.'
            ),
        ),
        'wave_delta': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'The wavelength step to use when generating the WCS (in Angstroms). '
                'If None, the default is set by the wavelength solution.'
            ),
        ),
        'astrometric': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If true, an astrometric correction will be applied using the alignment frames.'
            ),
        ),
        'scale_corr': newparset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=(
                'This option performs a small correction for the relative spectral illumination '
                'scale of different spec2D files. Specify the relative path+file to the spec2D '
                'file that you would like to use for the relative scaling. If you want to perform '
                'this correction, it is best to use the spec2d file with the highest S/N sky spectrum. '
                'You should choose the same frame for both the standards and science frames.'
            ),
        ),
        'correct_dar': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If True, the data will be corrected for differential atmospheric refraction (DAR).'
            ),
        ),
        'skysub_frame': newparset.set_parameter_definition(
            dtype=str,
            default='image',
            descr=(
                'Set the sky subtraction to be implemented. The default behaviour is to subtract '
                'the sky using the model that is derived from each individual image (i.e. set '
                'this parameter to "image"). To turn off sky subtraction completely, set this '
                'parameter to "none" (all lowercase). Finally, if you want to use a different frame '
                'for the sky subtraction, specify the relative path+file to the spec2D file that you '
                'would like to use for the sky subtraction. The model fit to the sky of the specified '
                'frame will be used. Note, the sky and science frames do not need to have the same '
                'exposure time; the sky model will be scaled to the science frame based on the relative exposure time.'
            ),
        ),
    }

    def validate(self):
        # Check the skysub options
        allowed_skysub_options = ["none", "image", ""]
        if self.data['skysub_frame'] not in allowed_skysub_options:
            if not Path(self.data['skysub_frame']).absolute().is_file():
                raise ValueError("The 'skysub_frame' must be one of:\n" + ", ".join(allowed_skysub_options) +
                                 "\nor, the relative path to a spec2d file.")
        if len(self.data['whitelight_range']) != 2:
            raise ValueError("The 'whitelight_range' must be a two element list of either NoneType or float")


class NewFluxCalibratePar(newparset.NewParSet):
    """
    New-style parameter set holding the arguments for how to perform the flux
    calibration (replacement for FluxCalibratePar).

    Mirrors the legacy `FluxCalibratePar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'fluxcalibrate'

    parameters = {
        'extrap_sens': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                "If False (default), the code will crash if one tries to use "
                "sensfunc at wavelengths outside its defined domain. By changing the "
                "par['sensfunc']['extrap_blu'] and par['sensfunc']['extrap_red'] this domain "
                "can be extended. If True the code will blindly extrapolate."
            ),
        ),
        'extinct_correct': newparset.set_parameter_definition(
            dtype=bool,
            default=None,
            descr=(
                'The default behavior for atmospheric extinction corrections is that if UVIS algorithm is used '
                '(which does not correct for telluric absorption) than an atmospheric extinction model '
                'is used to correct for extinction below 10,000A, whereas if the IR algorithm is used, then '
                'no extinction correction is applied since the atmosphere is modeled directly. To follow these '
                'defaults based on the algorithm this parameter should be set to ``extinct_correct=None``. If instead this '
                'parameter is set, this overide this default behavior. In other words, it will force an extinction correction '
                'if ``extinct_correct=True``, and will not perform an extinction correction if ``extinct_correct=False``.'
            ),
        ),
        'extinct_file': newparset.set_parameter_definition(
            dtype=str,
            default='closest',
            descr=(
                'If ``extinct_file=\'closest\'`` the code will select the PypeIt-included extinction '
                'file for the closest observatory (within 5 deg, geographic coordinates) to the telescope '
                'identified in ``std_file`` (see :ref:`extinction_correction` for the list of currently '
                'included files).  If constructing a sesitivity function for a telescope not within 5 deg '
                'of a listed observatory, this parameter may be set to the name of one of the listed '
                'extinction files.  Alternatively, a custom extinction file may be installed in the '
                'PypeIt cache using the ``pypeit_install_extinctfile`` script; this parameter may then '
                'be set to the name of the custom extinction file.'
            ),
        ),
        'use_archived_sens': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Use an archived sensfunc to flux calibration',
        ),
    }


class NewSensfuncUVISPar(newparset.NewParSet):
    """
    New-style parameter set for sensitivity function computation using the UV algorithm
    (replacement for SensfuncUVISPar).

    Mirrors the legacy `SensfuncUVISPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'sensfunc_uvis'

    parameters = {
        'std_file': newparset.set_parameter_definition(
            dtype=str,
            default=None,
            descr='Standard star file to generate sensfunc',
        ),
        'std_obj_id': newparset.set_parameter_definition(
            dtype=[str, int],
            default=None,
            descr=('Specifies object in spec1d file to use as standard. The brightest object found is used otherwise.'),
        ),
        'sensfunc': newparset.set_parameter_definition(
            dtype=str,
            default=None,
            descr='FITS file that contains or will contain the sensitivity function.',
        ),
        'extinct_correct': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If ``extinct_correct=True`` the code will use an atmospheric extinction model to '
                'extinction correct the data below 10000A. Note that this correction makes no '
                'sense if one is telluric correcting and this shold be set to False'
            ),
        ),
        'extinct_file': newparset.set_parameter_definition(
            dtype=str,
            default='closest',
            descr=(
                'If ``extinct_file=\'closest\'`` the code will select the PypeIt-included extinction '
                'file for the closest observatory (within 5 deg, geographic coordinates) to the telescope '
                'identified in ``std_file`` (see :ref:`extinction_correction` for the list of currently '
                'included files).  If constructing a sesitivity function for a telescope not within 5 deg '
                'of a listed observatory, this parameter may be set to the name of one of the listed '
                'extinction files.  Alternatively, a custom extinction file may be installed in the '
                'PypeIt cache using the ``pypeit_install_extinctfile`` script; this parameter may then '
                'be set to the name of the custom extinction file.'
            ),
        ),
        'telluric_correct': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                "If ``telluric_correct=True`` the code will grab the sens_dict['telluric'] tag from the "
                "sensfunc dictionary and apply it to the data."
            ),
        ),
        'telluric': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'If ``telluric=True`` the code creates a synthetic standard star spectrum using the Kurucz models, '
                'the sens func is created setting nresln=1.5 it contains the correction for telluric lines.'
            ),
        ),
        'polycorrect': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Whether you want to correct the sensfunc with polynomial in the telluric and recombination line regions',
        ),
        'polyfunc': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Whether you want to use the polynomial fit as your final SENSFUNC',
        ),
        'nresln': newparset.set_parameter_definition(
            dtype=[int, float],
            default=20,
            descr='Parameter governing the spacing of the bspline breakpoints in terms of number of resolution elements.',
        ),
        'resolution': newparset.set_parameter_definition(
            dtype=[int, float],
            default=3000.0,
            descr='Expected resolution of the standard star spectrum. This should be measured from the data.',
        ),
        'trans_thresh': newparset.set_parameter_definition(
            dtype=float,
            default=0.9,
            descr=(
                'Parameter for selecting telluric regions which are masked. Locations below this '
                'transmission value are masked. If you have significant telluric absorption you should '
                'be using telluric.sensnfunc_telluric'
            ),
        ),
    }

    def validate(self):
        if self.data['sensfunc'] is not None and self.data['std_file'] is None and not Path(self.data['sensfunc']).absolute().is_file():
            raise ValueError('Provided sensitivity function does not exist: {0}.'.format(self.data['sensfunc']))


class NewTelluricPar(newparset.NewParSet):
    """
    New-style parameter set holding telluric-correction arguments (replacement for TelluricPar).

    Mirrors the legacy `TelluricPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'telluric'

    valid_teltype = ['pca', 'grid']

    parameters = {
        'telgridfile': newparset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=(
                'File with the telluric model spectra to use.  Generally, these do '
                'not need to be set; reasonable defaults are provided for each '
                'spectrograph.  Due to their size, the files are not included with '
                'the released pypeit package; instead the code downloads each file '
                'into your cache as needed.  If this parameter is set in your pypeit '
                'file, it can be the path to a local file (which must have the '
                'correct format), or it can be the name of the specific cache file to '
                'use (e.g., TellPCA_3000_26000_R10000.fits).'
            ),
        ),
        'tell_npca': newparset.set_parameter_definition(
            dtype=int,
            default=5,
            descr='Number of telluric PCA components used. Can be set to any number from 1 to 10.',
        ),
        'teltype': newparset.set_parameter_definition(
            dtype=str,
            default='pca',
            options=valid_teltype,
            descr=(
                'Method used to evaluate telluric models.  Options are ``pca`` or '
                '``grid``. The ``grid`` option uses a fixed grid of pre-computed '
                'HITRAN+LBLRTM atmospheric transmission models for each observatory, '
                'whereas the ``pca`` option uses principal components of a larger '
                'model grid to compute an accurate pseudo-telluric model with a much '
                'lighter telgridfile.'
            ),
        ),
        'sn_clip': newparset.set_parameter_definition(
            dtype=[int, float],
            default=30.0,
            descr=(
                'This adds an error floor to the variance, preventing too much '
                'rejection at high-S/N (i.e., standard stars, bright objects), using '
                'the function :func:`~pypeit.utils.clip_ivar`. A small erorr is added '
                'to the input variance so that the output variance will never give '
                'S/N greater than ``sn_clip``. This prevents overly aggressive '
                'rejection in high S/N ratio spectra that neverthless differ at a '
                'level greater than the formal S/N due to the fact that our telluric '
                'models are only good to about 3%.'
            ),
        ),
        'resln_guess': newparset.set_parameter_definition(
            dtype=[int, float],
            default=None,
            descr=(
                'A guess for the resolution of your spectrum expressed as '
                'lambda/dlambda. The resolution is fit explicitly as part of the '
                'telluric model fitting, but this guess helps determine the bounds '
                'for the optimization (see ``resln_frac_bounds``). If not provided, '
                'the wavelength sampling of your spectrum will be used and the '
                'resolution calculated using a typical sampling of 3 spectral pixels '
                'per resolution element.'
            ),
        ),
        'resln_frac_bounds': newparset.set_parameter_definition(
            dtype=tuple,
            default=(0.6, 1.4),
            descr=(
                'Bounds for the resolution fit optimization which is part of the '
                'telluric model.  This range is in units of ``resln_guess``, so the '
                'default would bound the spectral resolution fit to be within the '
                'range ``bounds_resln = (0.6*resln_guess, 1.4*resln_guess)``.'
            ),
        ),
        'pix_shift_bounds': newparset.set_parameter_definition(
            dtype=tuple,
            default=(-5.0, 5.0),
            descr='Bounds for the pixel shift optimization in the telluric model fit in units of pixels.  The atmosphere will be allowed to shift within this range during the fit.',
        ),
        'delta_coeff_bounds': newparset.set_parameter_definition(
            dtype=tuple,
            default=(-20.0, 20.0),
            descr='Parameters setting the polynomial coefficient bounds for sensfunc optimization.',
        ),
        'minmax_coeff_bounds': newparset.set_parameter_definition(
            dtype=tuple,
            default=(-5.0, 5.0),
            descr=(
                "Parameters setting the polynomial coefficient bounds for sensfunc "
                "optimization.  Bounds are currently determined as follows.  We "
                "compute an initial fit to the sensfunc in the "
                ":func:`~pypeit.core.telluric.init_sensfunc_model` function. That "
                "determines a set of coefficients. The bounds are then determined "
                "according to: "
                "``[(np.fmin(np.abs(this_coeff)*obj_params['delta_coeff_bounds'][0], "
                "obj_params['minmax_coeff_bounds'][0]), "
                "np.fmax(np.abs(this_coeff)*obj_params['delta_coeff_bounds'][1], "
                "obj_params['minmax_coeff_bounds'][1]))]``."
            ),
        ),
        'maxiter': newparset.set_parameter_definition(
            dtype=int,
            default=2,
            descr='Maximum number of iterations for the telluric + object model fitting.  The code performs multiple iterations rejecting outliers at each step.  The fit is then performed anew to the remaining good pixels.  For this reason if you run with the ``disp=True`` option, you will see that the f(x) loss function gets progressively better during the iterations.',
        ),
        'sticky': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Sticky parameter for the :func:`~pypeit.utils.djs_reject` algorithm '
                'for iterative model fit rejection.  If set to True then points '
                'rejected from a previous iteration are kept rejected, in other words '
                'the bad pixel mask is the OR of all previous iterations and rejected '
                'pixels accumulate.  If set to False, the bad pixel mask is the mask '
                'from the previous iteration, and if the model fit changes between '
                'iterations, points can alternate from being rejected to not '
                'rejected.  At present this code only performs optimizations with '
                'differential evolution and experience shows that sticky needs to be '
                'True in order for these to converge.  This is because the outliers '
                'can be so large that they dominate the loss function, and one never '
                'iteratively converges to a good model fit.  In other words, the '
                'deformations in the model between iterations with ``sticky=False`` '
                'are too small to approach a reasonable fit.'
            ),
        ),
        'lower': newparset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr=(
                'Lower rejection threshold in units of ``sigma_corr*sigma``, where '
                '``sigma`` is the formal noise of the spectrum, and sigma_corr is an '
                'empirically determined correction to the formal error. The '
                'distribution of input chi (defined by ``chi = (data - '
                'model)/sigma``) values is analyzed, and a correction factor to the '
                'formal error ``sigma_corr`` is returned which is multiplied into the '
                'formal errors. In this way, a rejection threshold of e.g. 3 sigma, '
                'will always correspond to roughly the same percentile.  This '
                'renormalization is performed with '
                ':func:`~pypeit.coadd1d.renormalize_errors` function, and guarantees '
                'that rejection is not too aggressive in cases where the empirical '
                'errors determined from the chi-distribution differ significantly '
                'from the formal noise which is used to determine chi.'
            ),
        ),
        'upper': newparset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr='Upper rejection threshold in units of ``sigma_corr*sigma``, where ``sigma`` is the formal noise of the spectrum, and ``sigma_corr`` is an empirically determined correction to the formal error. See ``lower`` for additional detail.',
        ),
        'seed': newparset.set_parameter_definition(
            dtype=int,
            default=777,
            descr='An initial seed for the differential evolution optimization, which is a random process.  The default is 777, which will be used to generate a unique seed for every order.  A specific seed is used because otherwise the random number generator will use the time for the seed, and the results will not be reproducible.',
        ),
        'tol': newparset.set_parameter_definition(
            dtype=float,
            default=1e-3,
            descr='Relative tolerance for converage of the differential evolution optimization. See `scipy.optimize.differential_evolution`_ for details.',
        ),
        'popsize': newparset.set_parameter_definition(
            dtype=int,
            default=30,
            descr='A multiplier for setting the total population size for the differential evolution optimization. See `scipy.optimize.differential_evolution`_ for details.',
        ),
        'recombination': newparset.set_parameter_definition(
            dtype=[int, float],
            default=0.7,
            descr='The recombination constant for the differential evolution optimization. This should be in the range between 0 and 1. See `scipy.optimize.differential_evolution`_ for details.',
        ),
        'polish': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='If True then differential evolution will perform an additional optimization at the end to polish the best fit at the end, which can improve the optimization slightly. See `scipy.optimize.differential_evolution`_ for details.',
        ),
        'disp': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Argument for `scipy.optimize.differential_evolution`_ that will display status messages to the screen indicating the status of the optimization.  See documentation for :class:`~pypeit.core.telluric.Telluric` for a description of the output and how to know if things are working well.',
        ),
        'only_orders': newparset.set_parameter_definition(
            dtype=[int, list, np.ndarray],
            default=None,
            descr='Order number, or list of order numbers if you only want to fit specific orders.',
        ),
        'objmodel': newparset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=('The object model to be used for telluric fitting. Currently the options are: ``qso``, ``star``, and ``poly``.  For ``qso``, you might need to set ``redshift`` and ``bal_wv_min_max``.  For ``star``, you must set ``star_type``, ``star_ra``, ``star_dec``, and ``star_mag``.  For ``poly``, you might need to set ``fit_wv_min_max`` and ``norder``.'),
        ),
        'redshift': newparset.set_parameter_definition(
            dtype=[int, float],
            default=0.0,
            descr='The redshift for the object model. This is currently only used by the QSO model.',
        ),
        'delta_redshift': newparset.set_parameter_definition(
            dtype=float,
            default=0.1,
            descr='Range within the redshift can be varied for telluric fitting, i.e. the code performs a bounded optimization within the redshift +- delta_redshift.',
        ),
        'pca_file': newparset.set_parameter_definition(
            dtype=str,
            default='qso_pca_1200_3100.fits',
            descr='Fits file containing quasar PCA model. Needed for the QSO model.  If you change the default, you might need to set ``pca_lower`` and ``pca_upper``.',
        ),
        'npca': newparset.set_parameter_definition(
            dtype=int,
            default=8,
            descr='Number of pca for the objmodel=qso qso PCA fit',
        ),
        'bal_wv_min_max': newparset.set_parameter_definition(
            dtype=[list, np.ndarray],
            default=None,
            descr='Min/max wavelength of broad absorption features. If there are several BAL features, the format for this mask is ``[wave_min_bal1, wave_max_bal1, wave_min_bal2, wave_max_bal2,...]``. These masked pixels will be ignored during the fitting.',
        ),
        'bounds_norm': newparset.set_parameter_definition(
            dtype=tuple,
            default=(0.1, 3.0),
            descr='Normalization bounds for scaling the initial object model.',
        ),
        'tell_norm_thresh': newparset.set_parameter_definition(
            dtype=[int, float],
            default=0.9,
            descr='Threshold of telluric absorption region',
        ),
        'pca_lower': newparset.set_parameter_definition(
            dtype=[int, float],
            default=1220.0,
            descr='Minimum wavelength for the qso pca model',
        ),
        'pca_upper': newparset.set_parameter_definition(
            dtype=[int, float],
            default=3100.0,
            descr='Maximum wavelength for the qso pca model',
        ),
        'mask_lyman_a': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Mask the blueward of Lyman-alpha line during the fitting?',
        ),
        'star_type': newparset.set_parameter_definition(
            dtype=str,
            default=None,
            descr='stellar type',
        ),
        'star_mag': newparset.set_parameter_definition(
            dtype=[float, int],
            default=None,
            descr='AB magnitude in V band',
        ),
        'star_ra': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr='Object right-ascension in decimal deg',
        ),
        'star_dec': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr='Object declination in decimal deg',
        ),
        'func': newparset.set_parameter_definition(
            dtype=str,
            default='legendre',
            descr='Polynomial model function',
        ),
        'model': newparset.set_parameter_definition(
            dtype=str,
            default='exp',
            descr='Types of polynomial model. Options are ``poly``, ``square``, ``exp`` corresponding to normal polynomial, squared polynomial, or exponentiated polynomial.',
        ),
        'polyorder': newparset.set_parameter_definition(
            dtype=int,
            default=3,
            descr='Order of the polynomial model fit',
        ),
        'fit_wv_min_max': newparset.set_parameter_definition(
            dtype=list,
            default=None,
            descr='Pixels within this mask will be used during the fitting. The format is the same with ``bal_wv_min_max``, but this mask is good pixel masks.',
        ),
    }

    def validate(self):
        if self.data['tell_npca'] < 1 or self.data['tell_npca'] > 10:
            raise ValueError('Invalid value {:d} for tell_npca '.format(self.data['tell_npca'])+
                             '(must be between 1 and 10).')

        self.data['teltype'] = self.data['teltype'].lower()
        if self.data['teltype'] not in self.valid_teltype:
            raise ValueError('Invalid teltype "{}"'.format(self.data['teltype'])+
                             ', valid options are: {}.'.format(', '.join(self.valid_teltype)))


class NewSensFuncPar(newparset.NewParSet):
    """
    New-style parameter set holding the arguments for sensitivity function computation
    using the UV algorithm (replacement for SensFuncPar).

    Mirrors the legacy `SensFuncPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'sensfunc'

    parameters = {
        'use_flat': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='If True, the flatfield spectrum will be used when computing the sensitivity function.',
        ),
        'extr': newparset.set_parameter_definition(
            dtype=str,
            default='OPT',
            descr=(
                'Extraction method to use for the sensitivity function.  Options are: '
                "'OPT' (optimal extraction), 'BOX' (boxcar extraction). Default is 'OPT'."
            ),
        ),
        'extrap_blu': newparset.set_parameter_definition(
            dtype=float,
            default=0.1,
            descr=(
                'Fraction of minimum wavelength coverage to grow the wavelength coverage of the '
                'sensitivitity function in the blue direction (`i.e.`, if the standard star spectrum '
                'cuts off at ``wave_min``) the sensfunc will be extrapolated to cover down to '
                ' (1.0 - ``extrap_blu``) * ``wave_min``'
            ),
        ),
        'extrap_red': newparset.set_parameter_definition(
            dtype=float,
            default=0.1,
            descr=(
                'Fraction of maximum wavelength coverage to grow the wavelength coverage of the '
                'sensitivitity function in the red direction (`i.e.`, if the standard star spectrum'
                'cuts off at ``wave_max``) the sensfunc will be extrapolated to cover up to '
                ' (1.0 + ``extrap_red``) * ``wave_max``'
            ),
        ),
        'samp_fact': newparset.set_parameter_definition(
            dtype=float,
            default=1.5,
            descr=(
                'Sampling factor to make the wavelength grid for sensitivity function finer or coarser. '
                'samp_fact > 1.0 oversamples (finer), samp_fact < 1.0 undersamples (coarser).'
            ),
        ),
        'multi_spec_det': newparset.set_parameter_definition(
            dtype=list,
            default=None,
            descr=(
                'List of detectors (identified by their string name, like ' \
                "'DET01') to splice together for multi-detector instruments "
                '(e.g. DEIMOS). It is assumed that there is *no* overlap in ' \
                'wavelength across detectors (might be ok if there is).  If ' \
                "entered as a list of integers, they should be converted to ' "
                'the detector name.  **Cannot be used with detector mosaics.**'
            ),
        ),
        'trim_std_pixs': newparset.set_parameter_definition(
            dtype=[list, tuple],
            default=None,
            descr=(
                'List or tuple of two integers specifying the number of pixels to trim'
                'from the start and end of the 1D standard star spectrum. '
                'Example: [10, 5] will trim 10 pixels from the start (blue side)'
                'and 5 pixels from the end (red side) of the spectrum. '
            ),
        ),
        'algorithm': newparset.set_parameter_definition(
            dtype=str,
            default='UVIS',
            options=['UVIS', 'IR'],
            descr=(
                "Specify the algorithm for computing the sensitivity function. The options are: "
                r" (1) UVIS = Should be used for data with :math:`\lambda < 7000` A. "
                "No detailed model of telluric absorption but corrects for atmospheric extinction. "
                r" (2) IR = Should be used for data with :math:`\lambda > 7000` A. "
                "Peforms joint fit for sensitivity function and telluric absorption using HITRAN models."
            ),
        ),
        'UVIS': newparset.set_parameter_definition(
            dtype=NewSensfuncUVISPar,
            default=NewSensfuncUVISPar(),
            descr='Parameters for the UVIS sensfunc algorithm',
        ),
        'IR': newparset.set_parameter_definition(
            dtype=NewTelluricPar,
            default=NewTelluricPar(),
            descr='Parameters for the IR sensfunc algorithm',
        ),
        'polyorder': newparset.set_parameter_definition(
            dtype=[int, list],
            default=5,
            descr='Polynomial order for sensitivity function fitting',
        ),
        'star_type': newparset.set_parameter_definition(
            dtype=str,
            default=None,
            descr='Spectral type of the standard star (for near-IR mainly)',
        ),
        'star_mag': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr='Magnitude of the standard star (for near-IR mainly)',
        ),
        'star_ra': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr='RA of the standard star. This will override values in the header (`i.e.`, if they are wrong or absent)',
        ),
        'star_dec': newparset.set_parameter_definition(
            dtype=float,
            default=None,
            descr='DEC of the standard star. This will override values in the header (`i.e.`, if they are wrong or absent)',
        ),
        'mask_hydrogen_lines': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Mask hydrogen Balmer, Paschen, Brackett, and Pfund recombination lines in the sensitivity function fit. '
                'A region equal to ``hydrogen_mask_wid`` on either side of the line center is masked.'
            ),
        ),
        'hydrogen_mask_wid': newparset.set_parameter_definition(
            dtype=float,
            default=10.0,
            descr='Mask width from line center for hydrogen recombination lines in Angstroms (total mask width is 2x this value).',
        ),
        'mask_helium_lines': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Mask certain ``HeII`` recombination lines prominent in O-type stars in the sensitivity function fit '
                'A region equal to 0.5 * ``hydrogen_mask_wid`` on either side of the line center is masked.'
            ),
        ),
    }

    def validate(self):
        # Validate extraction choice
        allowed_extractions = ['BOX', 'OPT']
        if self.data['extr'] not in allowed_extractions:
            raise ValueError("'extr' must be one of:\n" + ", ".join(allowed_extractions))

        # check trim_std_pixs format
        if self.data['trim_std_pixs'] is not None:
            if not isinstance(self.data['trim_std_pixs'], (list, tuple)) or len(self.data['trim_std_pixs']) != 2:
                raise ValueError("`trim_std_pixs` must be a list or tuple of two integers.")


class NewSlitMaskPar(newparset.NewParSet):
    """
    New-style parameter set holding the arguments for slitmask ingestion and object assignment

    Mirrors the legacy `SlitMaskPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'slitmask'

    parameters = {
        'obj_toler': newparset.set_parameter_definition(
            dtype=[int, float],
            default=1.0,
            descr=(
                'If slitmask design information is provided, and slit matching is performed '
                '(``use_maskdesign = True`` in ``EdgeTracePar``), this parameter provides '
                'the desired tolerance (arcsec) to match sources to targeted objects'
            ),
        ),
        'assign_obj': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='If SlitMask object was generated, assign RA,DEC,name to detected objects',
        ),
        'use_alignbox': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Use stars in alignment boxes to compute the slitmask offset. '
                'If this is set to ``True`` PypeIt will NOT compute '
                'the offset using ``snr_thrshd`` or ``bright_maskdef_id``'
            ),
        ),
        'snr_thrshd': newparset.set_parameter_definition(
            dtype=[int, float],
            default=50.0,
            descr=(
                'Objects detected above this S/N threshold will '
                'be used to compute the slitmask offset. This is the default behaviour for DEIMOS '
                ' unless ``slitmask_offset``, ``bright_maskdef_id`` or ``use_alignbox`` is set.'
            ),
        ),
        'slitmask_offset': newparset.set_parameter_definition(
            dtype=[int, float],
            default=None,
            descr=(
                'User-provided slitmask offset (pixels) from the position expected by '
                'the slitmask design. This is optional, and if set PypeIt will NOT compute '
                'the offset using ``snr_thrshd`` or ``bright_maskdef_id``.'
            ),
        ),
        'use_dither_offset': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Use the dither offset recorded in the header of science frames as the value '
                'of the slitmask offset. This is currently only available for Keck MOSFIRE '
                'reduction and it is set as the default for this instrument. If set PypeIt will '
                'NOT compute the offset using ``snr_thrshd`` or ``bright_maskdef_id``. '
                'However, it is ignored if ``slitmask_offset`` is provided. '
            ),
        ),
        'bright_maskdef_id': newparset.set_parameter_definition(
            dtype=int,
            default=None,
            descr=(
                '`maskdef_id` (corresponding e.g., to `dSlitId` and `Slit_Number` '
                'in the DEIMOS/LRIS and MOSFIRE slitmask design, respectively) of a '
                'slit containing a bright object that will be used to compute the '
                'slitmask offset. This parameter is optional and is ignored '
                'if ``slitmask_offset`` is provided.'
            ),
        ),
        'extract_missing_objs': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Force extraction of undetected objects at the location expected '
                  'from the slitmask design.',
        ),
        'missing_objs_fwhm': newparset.set_parameter_definition(
            dtype=[int, float],
            default=None,
            descr=(
                'Indicates the FWHM in arcsec for the force extraction of undetected objects. '
                'PypeIt will try to determine the FWHM from the flux profile '
                '(by using ``missing_objs_fwhm`` as initial guess). '
                'If the FWHM cannot be determined, ``missing_objs_fwhm`` will be assumed. '
                'If you do not want PypeIt to try to determine the FWHM set the '
                'parameter ``use_user_fwhm`` in ``ExtractionPar`` to True. '
                'If ``missing_objs_fwhm`` is ``None`` (which is the default) PypeIt will use '
                'the median FWHM of all the detected objects.'
            ),
        ),
        'missing_objs_boxcar_rad': newparset.set_parameter_definition(
            dtype=[int, float],
            default=1.0,
            descr='Indicates the boxcar radius in arcsec for the force '
                  'extraction of undetected objects. ',
        ),
    }


class NewReduxPar(newparset.NewParSet):
    """
    New-style parameter set for global reduction settings (replacement for ReduxPar).

    Mirrors the legacy `ReduxPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'redux'

    parameters = {
        'spectrograph': newparset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=(
                'Spectrograph that provided the data to be reduced.  '
                'See :ref:`instruments` for valid options.'
            ),
        ),
        'quicklook': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Run a quick look reduction? This is usually good if you want to quickly '
                'reduce the data (usually at the telescope in real time) to get an initial '
                'estimate of the data quality.'
            ),
        ),
        'detnum': newparset.set_parameter_definition(
            dtype=[int, list],
            default=None,
            descr=(
                'Restrict reduction to a list of detector indices. '
                'In case of mosaic reduction (currently only available for '
                'Gemini/GMOS and Keck/DEIMOS) ``detnum`` should be a list of '
                'tuples of the detector indices that are mosaiced together. '
                'E.g., for Gemini/GMOS ``detnum`` would be ``[(1,2,3)]`` and for '
                'Keck/DEIMOS it would be ``[(1, 5), (2, 6), (3, 7), (4, 8)]``'
            ),
        ),
        'slitspatnum': newparset.set_parameter_definition(
            dtype=[str, list],
            default=None,
            descr=(
                'Restrict reduction to a set of slit DET:SPAT values (closest slit is used). '
                'Example syntax -- slitspatnum = DET01:175,DET01:205 or MSC02:2234  If you are re-running the code, '
                '(i.e. modifying one slit) you *must* have the precise SPAT_ID index.'
            ),
        ),
        'maskIDs': newparset.set_parameter_definition(
            dtype=[str, int, list],
            default=None,
            descr=(
                'Restrict reduction to a set of slitmask IDs '
                'Example syntax -- ``maskIDs = 818006,818015`` '
                'This must be used with detnum (for now).'
            ),
        ),
        'sortroot': newparset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=(
                'A filename given to output the details of the sorted files.  If '
                'None, the default is the root name of the pypeit file.  If off, '
                'no output is produced.'
            ),
        ),
        'calwin': newparset.set_parameter_definition(
            dtype=[int, float],
            default=0,
            descr=(
                'The window of time in hours to search for calibration frames for a '
                'science frame'
            ),
        ),
        'ignore_bad_headers': newparset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Ignore bad headers (NOT recommended unless you know it is safe).',
        ),
        'scidir': newparset.set_parameter_definition(
            dtype=str,
            default='Science',
            descr='Directory relative to calling directory to write science files.',
        ),
        'qadir': newparset.set_parameter_definition(
            dtype=str,
            default='QA',
            descr='Directory relative to calling directory to write quality '
                  'assessment files.',
        ),
        'redux_path': newparset.set_parameter_definition(
            dtype=str,
            default=os.getcwd(),
            descr=(
                'Path to folder for performing reductions.  Default is the '
                'current working directory.'
            ),
        ),
        'chk_version': newparset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If True enforce strict PypeIt version checking to ensure that '
                'all files were created with the current version of PypeIt.  If '
                'set to False, the code will attempt to read out-of-date files '
                'and keep going.  Beware (!!) that this can lead to unforeseen '
                'bugs that either cause the code to crash or lead to erroneous '
                'results. I.e., you really need to know what you are doing if '
                'you set this to False!'
            ),
        ),
    }

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

