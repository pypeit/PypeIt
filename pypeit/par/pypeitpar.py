"""
Defines parameter sets used to set the behavior for core pypeit
functionality.

.. include:: ../include/links.rst
"""
from pathlib import Path

from configobj import ConfigObj
from IPython import embed
import numpy as np
import os

from pypeit import dataPaths
from pypeit import log
from pypeit import PypeItError
from pypeit.core import parse
from pypeit.core.framematch import FrameTypeBitMask
from pypeit.par import parset
from pypeit.par import util


class TelescopePar(parset.ParSet):
    """
    New-style parameter set for the salient properties of a telescope.

    Mirrors the legacy `TelescopePar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'telescope'

    valid_telescopes = [
        'AAT', 'GEMINI-N','GEMINI-S', 'KECK', 'SHANE', 'WHT', 'APF', 'TNG', 'VLT',
        'MAGELLAN', 'LBT', 'MMT', 'KPNO', 'NOT', 'P200', 'BOK', 'GTC', 'SOAR', 'NTT',
        'LDT', 'JWST', 'HILTNER', 'SUBARU'
    ]

    parameters = {
        'name': parset.set_parameter_definition(
            dtype=str,
            default='KECK',
            options=valid_telescopes,
            descr=(
                'Name of the telescope used to obtain the observations.  Options are: '
                f'{valid_telescopes}'
            ),
        ),
        'longitude': parset.set_parameter_definition(
            dtype=[int, float],
            descr='Longitude of the telescope on Earth in degrees.',
        ),
        'latitude': parset.set_parameter_definition(
            dtype=[int, float],
            descr='Latitude of the telescope on Earth in degrees.',
        ),
        'elevation': parset.set_parameter_definition(
            dtype=[int, float],
            descr='Elevation of the telescope in m',
        ),
        'fratio': parset.set_parameter_definition(
            dtype=[int, float],
            descr='f-ratio of the telescope',
        ),
        'diameter': parset.set_parameter_definition(
            dtype=[int, float],
            descr='Diameter of the telescope in m',
        ),
        'eff_aperture': parset.set_parameter_definition(
            dtype=[int, float],
            descr='Effective aperture of the telescope in m^2',
        ),
    }

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
        return (
            None if self['fratio'] is None or self['diameter'] is None
            else 206265/self['fratio']/self['diameter']/1e3
        )

    # TODO This method is a place holder until we can get effective apertures
    # for all of the telescopes. I did my best but could not find them all
    # online.
    def eff_aperture(self):
        """
        Return the effective aperture of the telescope in square meters.
        """
        return (
            np.pi*self['diameter']**2/4.0 if self['eff_aperture'] is None else self['eff_aperture']
        )



class ScatteredLightPar(parset.ParSet):
    """
    The parameter set used to hold arguments for modeling the scattered light.

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
        'method': parset.set_parameter_definition(
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
        'finecorr_method': parset.set_parameter_definition(
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
        'finecorr_pad': parset.set_parameter_definition(
            dtype=int,
            default=4,
            descr=(
                'Number of unbinned pixels by which to extend the slit edges by when masking the '
                'slits for the fine correction to the scattered light.'
            )
        ),
        'finecorr_order': parset.set_parameter_definition(
            dtype=int,
            default=2,
            descr=(
                'Polynomial order to use for the fine correction to the scattered light '
                'subtraction. It should be a low value.'
            )
        ),
        'finecorr_mask': parset.set_parameter_definition(
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


class ProcessImagesPar(parset.ParSet):
    """
    New-style parameter set for basic image processing using `parset.ParSet`.

    This replaces the old instance-driven __init__ with a class-level
    `parameters` specification. The `parset.ParSet` base class handles defaulting,
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
        'trim': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Trim the image to the detector supplied region',
        ),
        'apply_gain': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Convert the ADUs to electrons using the detector gain',
        ),
        'orient': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Orient the raw image into the PypeIt frame',
        ),
        'use_biasimage': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Use a bias image.  If True, one or more must be supplied in the PypeIt file.',
        ),
        'use_overscan': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Subtract off the overscan.  Detector *must* have one or code will crash.',
        ),
        'overscan_method': parset.set_parameter_definition(
            dtype=str,
            default='savgol',
            options=valid_overscan_methods,
            descr=(
                'Method used to fit the overscan. '
                f'Options are: {", ".join(valid_overscan_methods)}  Note: Method "polynomial" '
                'is identical to "chebyshev"; the former is deprecated and will be removed.'
            ),
        ),
        'overscan_par': parset.set_parameter_definition(
            dtype=[int, list],
            default=[5, 65],
            descr=(
                'Parameters for the overscan subtraction.  For '
                "'chebyshev' or 'polynomial', set overcan_par = order; "
                "for 'savgol', set overscan_par = order, window size ; "
                "for 'median', set overscan_par = None or omit the keyword."
            ),
        ),
        'correct_nonlinear': parset.set_parameter_definition(
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
        'use_darkimage': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Subtract off a dark image.  If True, one or more darks must be provided.',
        ),
        'dark_expscale': parset.set_parameter_definition(
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
        'use_pattern': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Subtract off a detector pattern. This pattern is assumed to be '
                'sinusoidal along one direction, with a frequency that is '
                'constant across the detector.'
            ),
        ),
        'subtract_continuum': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Subtract off the continuum level from an image. This parameter should only '
                'be set to True to combine arcs with multiple different lamps. '
                'For all other cases, this parameter should probably be False.'
            ),
        ),
        'subtract_scattlight': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Subtract off the scattered light from an image. This parameter should only '
                'be set to True for spectrographs that have dedicated methods to subtract '
                'scattered light. For all other cases, this parameter should be False.'
            ),
        ),
        'scattlight': parset.set_parameter_definition(
            dtype=ScatteredLightPar,
            default=ScatteredLightPar(),
            descr='Scattered light subtraction parameters.',
        ),
        'empirical_rn': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'If True, use the standard deviation in the overscan region to '
                'measure an empirical readnoise to use in the noise model.'
            ),
        ),
        'shot_noise': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Use the bias- and dark-subtracted image to calculate and include '
                'electron count shot noise in the image processing error budget'
            ),
        ),
        'noise_floor': parset.set_parameter_definition(
            dtype=float,
            default=0.0,
            descr=(
                'Impose a noise floor by adding the provided fraction of the '
                'bias- and dark-subtracted electron counts to the error budget.  '
                'E.g., a value of 0.01 means that the S/N of the counts in the '
                'image will never be greater than 100.'
            ),
        ),
        'use_pixelflat': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Use the pixel flat to make pixel-level corrections.  A '
                'pixelflat image must be provied.'
            ),
        ),
        'use_illumflat': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Use the illumination flat to correct for the illumination profile of each slit.',
        ),
        'use_specillum': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Use the relative spectral illumination profiles to correct '
                'the spectral illumination profile of each slit. This is '
                'primarily used for slicer IFUs.  To use this, you must set '
                '``slit_illum_relative=True`` in the ``flatfield`` parameter set!'
            ),
        ),
        'skip_write_2d': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Skip writing the 2D spectrum for science frames.  WARNING: '
                'This option should only be considered for reducing the volume '
                'of output data when processing large numbers of frames and only '
                'after ensuring the quality of the resulting reductions.'
            ),
        ),
        'spat_flexure_correct': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Correct slits, illumination flat, etc. for flexure',
        ),
        'spat_flexure_maxlag': parset.set_parameter_definition(
            dtype=int,
            default=20,
            descr='Maximum of possible spatial flexure correction, in pixels',
        ),
        'spat_flexure_sigdetect': parset.set_parameter_definition(
            dtype=[int, float],
            default=5.0,
            descr=(
                'Sigma threshold above fluctuations in the '
                'Sobel-filtered significance image, used for '
                'finding slit edges in the spectral image, '
                'for which the spatial flexure is computed.'
            ),
        ),
        'spat_flexure_vrange': parset.set_parameter_definition(
            dtype=tuple,
            descr=(
                'This parameter is used when generating the QA plot for the spatial flexure. '
                'It sets the data range (vmin,vmax) used by the colormap when showing the '
                'spectral image. If None, the range is set automatically.'
            ),
        ),
        'combine': parset.set_parameter_definition(
            dtype=str,
            default='mean',
            options=valid_combine_methods,
            descr=(
                'Method used to combine multiple frames.  Options are: '
                f'{", ".join(valid_combine_methods)}'
            ),
        ),
        'clip': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Perform sigma clipping when combining.  Only used with combine=mean',
        ),
        'scale_to_mean': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='If True, scale the input images to have the same mean before combining.',
        ),
        'comb_sigrej': parset.set_parameter_definition(
            dtype=float,
            descr=(
                'Sigma-clipping level for when clip=True; '
                'Use None for automatic limit (recommended).  '
            ),
        ),
        'satpix': parset.set_parameter_definition(
            dtype=str,
            default='reject',
            options=valid_saturation_handling,
            descr=(
                'Handling of saturated pixels.  Options are: '
                f'{", ".join(valid_saturation_handling)}'
            ),
        ),
        'mask_cr': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Identify CRs and mask them',
        ),
        'n_lohi': parset.set_parameter_definition(
            dtype=list,
            default=[0, 0],
            descr=(
                'Number of pixels to reject at the lowest and highest ends of the '
                'distribution; i.e., n_lohi = low, high.  Use None for no limit.'
            ),
        ),
        'lamaxiter': parset.set_parameter_definition(
            dtype=int,
            default=1,
            descr='Maximum number of iterations for LA cosmics routine.',
        ),
        'grow': parset.set_parameter_definition(
            dtype=[int, float],
            default=1.5,
            descr=(
                'Factor by which to expand regions with cosmic rays detected by the '
                'LA cosmics routine.'
            ),
        ),
        'rmcompact': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Remove compact detections in LA cosmics routine',
        ),
        'sigclip': parset.set_parameter_definition(
            dtype=[int, float],
            default=4.5,
            descr='Sigma level for rejection in LA cosmics routine',
        ),
        'sigfrac': parset.set_parameter_definition(
            dtype=[int, float],
            default=0.3,
            descr='Fraction for the lower clipping threshold in LA cosmics routine.',
        ),
        'objlim': parset.set_parameter_definition(
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


class FrameGroupPar(parset.ParSet):
    """
    The abstract base class for each frame type and the details of how they
    should be processed.
    """

    frametype = None
    """
    The frametype for this ParSet, which must be overwritten by the subclass
    """

    parameters = {
        'exprng': parset.set_parameter_definition(
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
    """
    The base-class only defines the exposure range for the frametype.
    """

    def validate(self):
        """
        Validate the frame parameters.

        The ``exprng`` must be valid, and the frametype must be one of the
        allowed values.
        """
        if self.data['exprng'] is not None and len(self.data['exprng']) != 2:
            raise ValueError('exprng must be a list with two items.')

        if self.frametype is None:
            raise ValueError('CODING ERROR: Subclasses of FrameGroupPar must define the frametype')
        valid_frametypes = FrameTypeBitMask().keys()
        if self.frametype not in valid_frametypes:
            raise ValueError(
                f'{self.frametype} is not a valid frametype.  Options are: {valid_frametypes}'
            )


class BiasFramePar(FrameGroupPar):
    frametype = 'bias'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
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


class DarkFramePar(FrameGroupPar):
    frametype = 'dark'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
                mask_cr=True,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class ScatteredLightFramePar(FrameGroupPar):
    frametype = 'scattlight'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                satpix='nothing',
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class PixelFlatFramePar(FrameGroupPar):
    frametype = 'pixelflat'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                satpix='nothing',
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class IllumFlatFramePar(FrameGroupPar):
    frametype = 'illumflat'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                satpix='nothing',
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class LampOffFlatsFramePar(FrameGroupPar):
    frametype = 'lampoffflats'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                satpix='nothing',
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class SlitlessPixFlatFramePar(FrameGroupPar):
    frametype = 'slitless_pixflat'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                satpix='nothing',
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
                combine='median',
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class PinholeFramePar(FrameGroupPar):
    frametype = 'pinhole'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(),
            descr='Low level parameters used for basic image processing',
        ),
    }


class AlignFramePar(FrameGroupPar):
    frametype = 'align'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                satpix='nothing',
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class ArcFramePar(FrameGroupPar):
    frametype = 'arc'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class TiltFramePar(FrameGroupPar):
    frametype = 'tilt'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class TraceFramePar(FrameGroupPar):
    frametype = 'trace'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                use_pixelflat=False,
                use_illumflat=False,
                use_specillum=False,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class StandardFramePar(FrameGroupPar):
    frametype = 'standard'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                noise_floor=0.01,
                mask_cr=True,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class SkyFramePar(FrameGroupPar):
    frametype = 'sky'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                noise_floor=0.01,
                mask_cr=True,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class ScienceFramePar(FrameGroupPar):
    frametype = 'science'
    default_key = f'{frametype}frame'

    parameters = FrameGroupPar.parameters | {
        'process': parset.set_parameter_definition(
            dtype=ProcessImagesPar,
            default=ProcessImagesPar(
                noise_floor=0.01,
                mask_cr=True,
            ),
            descr='Low level parameters used for basic image processing',
        ),
    }


class FlatFieldPar(parset.ParSet):
    """
    New-style parameter set for flat-fielding (replacement for FlatFieldPar).

    Mirrors the legacy `FlatFieldPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'flatfield'

    valid_methods = ['bspline', 'skip']

    valid_tweak_methods = ['threshold', 'gradient']
    
    valid_saturated_slits_methods = ['crash', 'mask', 'continue']

    parameters = {
        'method': parset.set_parameter_definition(
            dtype=str,
            default='bspline',
            options=valid_methods,
            descr=(
                'Method used to flat field the data; use skip to skip flat-fielding.  '
                f'Options are: None, {", ".join(valid_methods)}'
            ),
        ),
        'pixelflat_file': parset.set_parameter_definition(
            dtype=str,
            default=None,
            descr='Filename of the image to use for pixel-level field flattening',
        ),
        'spec_samp_fine': parset.set_parameter_definition(
            dtype=[int, float],
            default=1.2,
            descr=(
                'bspline break point spacing in units of pixels for spectral fit to '
                'flat field blaze function.'
            ),
        ),
        'spec_samp_coarse': parset.set_parameter_definition(
            dtype=[int, float],
            default=50.0,
            descr=(
                'bspline break point spacing in units of pixels for 2-d '
                'bspline-polynomial fit to flat field image residuals. '
                'This should be a large number unless you are trying to '
                'fit a sky flat with lots of narrow spectral features.'
            ),
        ),
        'spat_samp': parset.set_parameter_definition(
            dtype=[int, float],
            default=5.0,
            descr=(
                'Spatial sampling for slit illumination function. This is the width of the '
                'median filter in pixels used to determine the slit illumination function, '
                'and thus sets the minimum scale on which the illumination function will '
                'have features.'
            ),
        ),
        'pixelflat_min_wave': parset.set_parameter_definition(
            dtype=[int, float],
            default=None,
            descr=(
                'All values of the normalized pixel flat are set to 1 for '
                'wavelengths below this value.'
            ),
        ),
        'pixelflat_max_wave': parset.set_parameter_definition(
            dtype=[int, float],
            default=None,
            descr=(
                'All values of the normalized pixel flat are set to 1 for '
                'wavelengths above this value.'
            ),
        ),
        'tweak_slits': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Use the illumination flat field to tweak the slit edges. '
                'This will work even if illumflatten is set to False '
            ),
        ),
        'tweak_method': parset.set_parameter_definition(
            dtype=str,
            default='threshold',
            options=valid_tweak_methods,
            descr=(
                'Method used to tweak the slit edges (when "tweak_slits" is set to True). '
                f'Options include: {", ".join(valid_tweak_methods)}. '
                'The "threshold" method determines when the left and right slit edges '
                'fall below a threshold relative to the peak illumination. '
                'The "gradient" method determines where the gradient is the highest at '
                'the left and right slit edges. This method performs better when there is '
                'systematic vignetting in the spatial direction.'
            ),
        ),
        'tweak_slits_thresh': parset.set_parameter_definition(
            dtype=float,
            default=0.93,
            descr=(
                'If tweak_slits is True, this sets the illumination function threshold used to '
                'tweak the slit boundaries based on the illumination flat. '
                'It should be a number less than 1.0'
            ),
        ),
        'tweak_slits_maxfrac': parset.set_parameter_definition(
            dtype=float,
            default=0.10,
            descr=(
                'If tweak_slit is True, this sets the maximum fractional amount '
                '(of a slits width) allowed for trimming each (i.e. left and right) '
                'slit boundary, i.e. the default is 10% which means slits would '
                'shrink or grow by at most 20% (10% on each side)'
            ),
        ),
        'rej_sticky': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Propagate the rejected pixels through the stages of the '
                'flat-field fitting (i.e, from the spectral fit, to the spatial '
                'fit, and finally to the 2D residual fit).  If False, pixels '
                'rejected in each stage are included in each subsequent stage.'
            ),
        ),
        'slit_trim': parset.set_parameter_definition(
            dtype=[int, float, tuple],
            default=3.0,
            descr=(
                'The number of pixels to trim each side of the slit when '
                'selecting pixels to use for fitting the spectral response '
                'function.  Single values are used for both slit edges; a '
                'two-tuple can be used to trim the left and right sides differently.'
            ),
        ),
        'slit_illum_pad': parset.set_parameter_definition(
            dtype=[int, float],
            default=5.0,
            descr=(
                'The number of pixels to pad the slit edges when constructing the '
                'slit-illumination profile. Single value applied to both edges.'
            ),
        ),
        'slit_illum_finecorr': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If True, a fine correction to the spatial illumination profile '
                'will be performed. The fine correction is a low order 2D polynomial '
                'fit to account for a gradual change to the spatial illumination '
                'profile as a function of wavelength.'
            ),
        ),
        'slit_illum_relative': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Generate an image of the relative spectral illumination '
                'for a multi-slit setup.  If you set ``use_specillum = '
                'True`` for any of the frames that use the flatfield '
                'model, this *must* be set to True. Currently, this is '
                'only used for SlicerIFU reductions.'
            ),
        ),
        'illum_iter': parset.set_parameter_definition(
            dtype=int,
            default=0,
            descr=(
                'The number of rejection iterations to perform when constructing the '
                'slit-illumination profile.  No rejection iterations are performed '
                'if 0.  WARNING: Functionality still being tested.'
            ),
        ),
        'illum_rej': parset.set_parameter_definition(
            dtype=[int, float],
            default=5.0,
            descr=(
                'The sigma threshold used in the rejection iterations used to refine '
                'the slit-illumination profile.  Rejection iterations are only '
                'performed if ``illum_iter > 0``.'
            ),
        ),
        'twod_fit_npoly': parset.set_parameter_definition(
            dtype=int,
            default=None,
            descr=(
                'Order of polynomial used in the 2D bspline-polynomial fit to '
                'flat-field image residuals. The code determines the order of '
                'these polynomials to each slit automatically depending on '
                'the slit width, which is why the default is None. Alter '
                'this paramter at your own risk!'
            ),
        ),
        'saturated_slits': parset.set_parameter_definition(
            dtype=str,
            default='crash',
            options=valid_saturated_slits_methods,
            descr=(
                'Behavior when a slit is encountered with a large fraction '
                'of saturated pixels in the flat-field.  The options are: '
                "'crash' - Raise an error and halt the data reduction; "
                "'mask' - Mask the slit, meaning no science data will be "
                "extracted from the slit; 'continue' - ignore the "
                'flat-field correction, but continue with the reduction.'
            ),
        ),
        'slit_illum_ref_idx': parset.set_parameter_definition(
            dtype=int,
            default=0,
            descr=(
                'The index of a reference slit (0-indexed) used for estimating the '
                'relative spectral sensitivity (or the relative blaze). This parameter '
                'is only used if ``slit_illum_relative = True``.'
            ),
        ),
        'slit_illum_smooth_npix': parset.set_parameter_definition(
            dtype=int,
            default=10,
            descr=(
                'The number of pixels used to determine smoothly varying relative '
                'weights is given by ``nspec/slit_illum_smooth_npix``, where nspec is '
                'the number of spectral pixels.'
            ),
        ),
        'fit_2d_det_response': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Set this variable to True if you want to compute and '
                'account for the detector response in the flatfield image. '
                'Note that ``detector response`` refers to pixel sensitivity '
                'variations that primarily depend on (x,y) detector coordinates. '
                'In most cases, the default 2D bspline is sufficient to account '
                'for detector response (i.e. set this parameter to False). Note '
                'that this correction will _only_ be performed for the spectrographs '
                'that have a dedicated response correction implemented. Currently,'
                'this correction is only implemented for Keck+KCWI.'
            ),
        ),
    }

    def validate(self):
        """
        Check the parameters are valid for the provided method.
        """
        if self.data['pixelflat_file'] is None:
            return

        # Check the frame exists
        file_path = dataPaths.pixelflat.get_file_path(
            self.data['pixelflat_file'], return_none=True
        )
        if file_path is None:
            raise PypeItError(
                f'Provided pixelflat file, {self.data["pixelflat_file"]} not found. It is not a '
                'direct path, a cached file, or a file that can be downloaded from a PypeIt '
                'repository.'
            )


class FlexurePar(parset.ParSet):
    """
    New-style parameter set for flexure correction parameters.

    Mirrors the legacy `FlexurePar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'flexure'

    valid_methods = ['boxcar', 'slitcen', 'skip']

    valid_excessive_shift_methods = ['crash', 'set_to_zero', 'continue', 'use_median']

    parameters = {
        'spec_method': parset.set_parameter_definition(
            dtype=str,
            default='skip',
            options=valid_methods,
            descr=(
                'Method used to correct for flexure. Use skip for no correction.  If '
                'slitcen is used, the flexure correction is performed before the '
                'extraction of objects (not recommended).  '
                f'Options are: {", ".join(valid_methods)}'
            ),
        ),
        'spec_maxshift': parset.set_parameter_definition(
            dtype=int,
            default=20,
            descr='Maximum allowed spectral flexure shift in pixels.',
        ),
        'spectrum': parset.set_parameter_definition(
            dtype=str,
            default='paranal_sky.fits',
            descr=(
                'Archive sky spectrum to be used for the flexure correction. '
                'See ``pypeit/data/sky_spec/`` for a list of available sky spectra. '
                'If ``model`` is used, a model sky spectrum will be generated '
                'using :func:`~pypeit.wavemodel.nearIR_modelsky` and the spectral'
                'resolution of the spectrum to be flexure corrected.'
            ),
        ),
        'excessive_shift': parset.set_parameter_definition(
            dtype=str,
            default='use_median',
            options=valid_excessive_shift_methods,
            descr=(
                'Behavior when the measured spectral flexure shift is '
                'larger than ``spec_maxshift``.  The options are: '
                "'crash' - Raise an error and halt the data reduction; "
                "'set_to_zero' - Set the flexure shift to zero and continue "
                "with the reduction; 'continue' - Use the large "
                "flexure value whilst issuing a warning; and 'use_median' - "
                "Use the median flexure shift among all the objects in the same slit "
                "(if more than one object is detected) or among all "
                "the other slits; if not available, the flexure correction will not be applied."
            ),
        ),
        'minwave': parset.set_parameter_definition(
            dtype=[int, float],
            default=None,
            descr=(
                'Minimum wavelength to use for the correlation.  If ``None`` or less than '
                'the minimum wavelength of either the object or archive sky spectrum, this '
                'this parameter has no effect.'
            ),
        ),
        'maxwave': parset.set_parameter_definition(
            dtype=[int, float],
            default=None,
            descr=(
                'Maximum wavelength to use for the correlation.  If ``None`` or greater than '
                'the maximum wavelength of either the object or archive sky spectrum, this '
                'this parameter has no effect.'
            ),
        ),
        'multi_min_SN': parset.set_parameter_definition(
            dtype=[int, float],
            default=1,
            descr='Minimum S/N for analyzing sky spectrum for flexure',
        ),
    }


class AlignPar(parset.ParSet):
    """
    New-style parameter set for alignment tracing (replacement for AlignPar).

    Mirrors the legacy `AlignPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'align'

    parameters = {
        'locations': parset.set_parameter_definition(
            dtype=[list, np.ndarray],
            default=[0.0, 1.0],
            descr='Locations of the bars, in a list, specified as a fraction of the slit width',
        ),
        'trace_npoly': parset.set_parameter_definition(
            dtype=int,
            default=4,
            descr='Order of the polynomial to use when fitting the trace of a single bar',
        ),
        'trim_edge': parset.set_parameter_definition(
            dtype=list,
            default=[0, 0],
            descr='Trim the slit by this number of pixels left/right before finding alignment bars',
        ),
        'snr_thresh': parset.set_parameter_definition(
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


class Coadd1DPar(parset.ParSet):
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
        'ex_value': parset.set_parameter_definition(
            dtype=str,
            default='OPT',
            options=valid_extractions,
            descr="The extraction to coadd, i.e. optimal or boxcar. Must be either 'OPT' or 'BOX'",
        ),
        'flux_value': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='If True (default), the code will coadd the fluxed spectra (i.e. the FLAM) in the spec1d files. If False, it will coadd the counts.',
        ),
        'nmaskedge': parset.set_parameter_definition(
            dtype=int,
            default=2,
            descr='Number of edge pixels to mask. This should be removed/fixed.',
        ),
        'sn_smooth_npix': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Number of pixels to median filter by when computing S/N used to decide how to scale '
                'and weight spectra. If set to None (default), the code will determine the effective '
                'number of good pixels per spectrum in the stack that is being co-added and use 10% of '
                'this neff.'
            ),
        ),
        'sigrej_exp': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Rejection threshold used for rejecting exposures with S/N more than sigrej_exp*sigma '
                'above the median S/N. If None (the default), no rejection is performed. Currently, '
                'only available for multi-slit observations.'
            ),
        ),
        'wave_method': parset.set_parameter_definition(
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
        'dv': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                "Dispersion in units of km/s in case you want to specify it in the get_wave_grid  (for the 'velocity' option), "
                "otherwise a median value is computed from the data."
            ),
        ),
        'dwave': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                "Dispersion in Angstroms in case you want to specify it in the get_wave_grid  (for the 'linear' option), "
                "otherwise a median value is computed from the data."
            ),
        ),
        'dloglam': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                "Dispersion in units of log10(wave) in case you want to specify it in the get_wave_grid  (for the 'velocity' or 'log10' options), "
                "otherwise a median value is computed from the data."
            ),
        ),
        'wave_grid_min': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Used in case you want to specify the minimum wavelength in your wavelength grid, default=None computes from data'
            ),
        ),
        'wave_grid_max': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Used in case you want to specify the maximum wavelength in your wavelength grid, default=None computes from data'
            ),
        ),
        'spec_samp_fact': parset.set_parameter_definition(
            dtype=float,
            default=1.0,
            descr=(
                "Make the wavelength grid  sampling finer (spec_samp_fact < 1.0) or coarser "
                "(spec_samp_fact > 1.0) by this sampling factor. This basically multiples the 'native' "
                "spectral pixels by spec_samp_fact, i.e. units spec_samp_fact are pixels."
            ),
        ),
        'ref_percentile': parset.set_parameter_definition(
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
        'maxiter_scale': parset.set_parameter_definition(
            dtype=int,
            default=5,
            descr='Maximum number of iterations performed for rescaling spectra.',
        ),
        'sigrej_scale': parset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr='Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec.',
        ),
        'scale_method': parset.set_parameter_definition(
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
        'sn_min_medscale': parset.set_parameter_definition(
            dtype=[int, float],
            default=0.5,
            descr='For scale method set to ``auto``, this sets the minimum SNR for which median scaling is attempted.',
        ),
        'sn_min_polyscale': parset.set_parameter_definition(
            dtype=[int, float],
            default=2.0,
            descr='For scale method set to ``auto``, this sets the minimum SNR for which polynomial scaling is attempted.',
        ),
        'weight_method': parset.set_parameter_definition(
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
        'maxiter_reject': parset.set_parameter_definition(
            dtype=int,
            default=5,
            descr=(
                'Maximum number of iterations for stacking and rejection. The code stops iterating '
                'either when the output mask does not change betweeen successive iterations or when '
                'maxiter_reject is reached.'
            )
        ),
        'lower': parset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr='Lower rejection threshold used for rejecting pixels when combining spectra in units of sigma.',
        ),
        'upper': parset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr='Upper rejection threshold used for rejecting pixels when combining spectra in units of sigma.',
        ),
        'maxrej': parset.set_parameter_definition(
            dtype=int,
            descr=(
                'Coadding performs iterative rejection by comparing each exposure to a preliminary stack of '
                'all the exposures. If this parameter is set then it will not reject more than maxrej pixels '
                'per iteration of this rejection. The default is None, which means no maximum on rejected pixels.'
            ),
        ),
        'sn_clip': parset.set_parameter_definition(
            dtype=[int, float],
            default=30.0,
            descr=(
                'Errors are capped during rejection so that the S/N is never greater than sn_clip. This '
                'prevents overly aggressive rejection in high S/N ratio spectrum which neverthless differ '
                'at a level greater than the formal S/N due to systematics.'
            ),
        ),
        'nbests': parset.set_parameter_definition(
            dtype=[list, int],
            descr=(
                'Number of orders to use for estimating the per exposure weights. Default is None, '
                'which will just use one fourth of the total number of orders. This is only used for Echelle'
            ),
        ),
        'filter': parset.set_parameter_definition(
            dtype=str,
            default='none',
            descr='Filter for scaling.  See flux_calib.load_fitler_file() for naming.  Ignore if none',
        ),
        'mag_type': parset.set_parameter_definition(
            dtype=str,
            default='AB',
            descr='Magnitude type.  AB is the only option currently allowed',
        ),
        'filter_mag': parset.set_parameter_definition(
            dtype=float,
            descr='Magnitude of the source in the given filter',
        ),
        'filter_mask': parset.set_parameter_definition(
            dtype=[str, list],
            descr=(
                'List of wavelength regions to mask when doing the scaling (`i.e.`, occasional junk pixels). '
                'Colon and comma separateed, e.g.   5552:5559,6010:6030'
            ),
        ),
        'coaddfile': parset.set_parameter_definition(
            dtype=str,
            descr='Output filename',
        ),
    }


class Coadd2DPar(parset.ParSet):
    """
    New-style parameter set for 2D coaddition (replacement for Coadd2DPar).

    Mirrors the legacy `Coadd2DPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'coadd2d'

    parameters = {
        'only_slits': parset.set_parameter_definition(
            dtype=[str, list],
            default=None,
            descr=(
                'Restrict coaddition to one or more of slits. Example syntax -- '
                'DET01:175,DET02:205 or MSC02:2234. This and ``exclude_slits`` '
                'are mutually exclusive. If both are provided, ``only_slits`` takes precedence.'
            ),
        ),
        'exclude_slits': parset.set_parameter_definition(
            dtype=[str, list],
            default=None,
            descr=(
                'Exclude one or more slits from the coaddition. Example syntax -- '
                'DET01:175,DET02:205 or MSC02:2234. This and ``only_slits`` '
                'are mutually exclusive. If both are provided, ``only_slits`` takes precedence.'
            ),
        ),
        'offsets': parset.set_parameter_definition(
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
        'spat_toler': parset.set_parameter_definition(
            dtype=int,
            default=5,
            descr=(
                'This parameter provides the desired tolerance in spatial pixel used to identify '
                'slits in different exposures'
            )
        ),
        'use_slits4wvgrid': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='If True, use the slits to set the trace down the center',
        ),
        'weights': parset.set_parameter_definition(
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
        'user_obj_ids': parset.set_parameter_definition(
            dtype=list,
            default=None,
            descr=(
                'List of unique object identifiers that the user wants to use '
                'to compute the weights and/or the offsets for coadding images. '
                'For longslit/multislit spectroscopy, provide the ``SPAT_PIXPOS_ID`` '
                'of the object in each of the exposures. For echelle spectroscopy, '
                'provide the ``ECH_FRACPOS_ID`` of the object in each exposure. '
                'These unique object identifiers can be found in the spec1d*.txt '
                'files for each exposure. See :doc:`out_spec1D` for more info about '
                '``SPAT_PIXPOS_ID`` and ``ECH_FRACPOS_ID``. This parameter must always '
                'be a list of the same length as the number of exposures being coadded. '
                'If this parameter is not ``None``, it will be used to compute the offsets '
                'only if ``offsets = auto``, and it will used to compute the weights '
                'only if ``weights = auto``.'
            ),
        ),
        'manual': parset.set_parameter_definition(
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
        'wave_method': parset.set_parameter_definition(
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
        'spec_samp_fact': parset.set_parameter_definition(
            dtype=float,
            default=1.0,
            descr=(
                "Make the wavelength grid sampling finer (``spec_samp_fact`` less than 1.0)"
                "or coarser (``spec_samp_fact`` greater than 1.0) by this sampling factor."
                "This  multiples the 'native' spectral pixel size by ``spec_samp_fact``,"
                "i.e. the units of ``spec_samp_fact`` are pixels."
            ),
        ),
        'spat_samp_fact': parset.set_parameter_definition(
            dtype=float,
            default=1.0,
            descr=(
                "Make the spatial sampling finer (``spat_samp_fact`` less"
                "than 1.0) or coarser (``spat_samp_fact`` greather than 1.0) by"
                "this sampling factor. This basically multiples the 'native'"
                "spatial pixel size by ``spat_samp_fact``, i.e. the units of"
                "``spat_samp_fact`` are pixels."
            ),
        ),
    }

    def validate(self):
        # Normalize manual extraction entries if provided
        if self.data['manual'] is not None:
            self.data['manual'] = ';'.join(parse.fix_config_par_image_location(self.data['manual']))


class CubePar(parset.ParSet):
    """
    New-style parameter set for cube generation (replacement for CubePar).

    Mirrors the legacy `CubePar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'cube'

    valid_weight_methods = ['auto', 'constant', 'uniform', 'wave_dependent', 'relative', 'ivar']

    parameters = {
        'slit_spec': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If the data use slits in one spatial direction, set this to True. '
                'If the data uses fibres for all spaxels, set this to False.'
            ),
        ),
        'weight_method': parset.set_parameter_definition(
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
        'align': parset.set_parameter_definition(
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
        'combine': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='If set to True, the input frames will be combined. Otherwise, a separate datacube will be generated for each input spec2d file, and will be saved as a spec3d file.',
        ),
        'output_filename': parset.set_parameter_definition(
            dtype=str,
            default="",
            descr=(
                'If combining multiple frames, this string sets the output filename of '
                'the combined datacube. If combine=False, the output filenames will be '
                'prefixed with ``spec3d_*``'
            ),
        ),
        'sensfile': parset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=(
                'Filename of a sensitivity function to use to flux calibrate your datacube. '
                'The sensitivity function file will also be used to correct the relative scales '
                'of the slits.'
            ),
        ),
        'reference_image': parset.set_parameter_definition(
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
        'save_whitelight': parset.set_parameter_definition(
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
        'whitelight_range': parset.set_parameter_definition(
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
        'method': parset.set_parameter_definition(
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
        'spec_subpixel': parset.set_parameter_definition(
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
        'spat_subpixel': parset.set_parameter_definition(
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
        'slice_subpixel': parset.set_parameter_definition(
            dtype=int,
            default=5,
            descr=(
                'When method=subpixel, slice_subpixel sets the subpixellation scale of '
                'each IFU slice. The default option is to divide each slice into 5 sub-slices '
                'during datacube creation. See also, spec_subpixel and spat_subpixel.'
            ),
        ),
        'ra_min': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Minimum RA to use when generating the WCS. If None, the default is minimum RA '
                'based on the WCS of all spaxels. Units should be degrees.'
            ),
        ),
        'ra_max': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Maximum RA to use when generating the WCS. If None, the default is maximum RA '
                'based on the WCS of all spaxels. Units should be degrees.'
            ),
        ),
        'dec_min': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Minimum DEC to use when generating the WCS. If None, the default is minimum DEC '
                'based on the WCS of all spaxels. Units should be degrees.'
            ),
        ),
        'dec_max': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Maximum DEC to use when generating the WCS. If None, the default is maximum DEC '
                'based on the WCS of all spaxels. Units should be degrees.'
            ),
        ),
        'wave_min': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Minimum wavelength to use when generating the WCS. If None, the default is '
                'minimum wavelength based on the WCS of all spaxels. Units should be Angstroms.'
            ),
        ),
        'wave_max': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'Maximum wavelength to use when generating the WCS. If None, the default is '
                'maximum wavelength based on the WCS of all spaxels. Units should be Angstroms.'
            ),
        ),
        'spatial_delta': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'The spatial size of each spaxel to use when generating the WCS (in arcsec). '
                'If None, the default is set by the spectrograph file.'
            ),
        ),
        'wave_delta': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'The wavelength step to use when generating the WCS (in Angstroms). '
                'If None, the default is set by the wavelength solution.'
            ),
        ),
        'astrometric': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If true, an astrometric correction will be applied using the alignment frames.'
            ),
        ),
        'scale_corr': parset.set_parameter_definition(
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
        'correct_dar': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If True, the data will be corrected for differential atmospheric refraction (DAR).'
            ),
        ),
        'skysub_frame': parset.set_parameter_definition(
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
                raise ValueError(
                    'The "skysub_frame" must be one of the following options: '
                    f'{", ".join(allowed_skysub_options)} or, the relative path to a spec2d file.'
                )
        if len(self.data['whitelight_range']) != 2:
            raise ValueError(
                "The 'whitelight_range' must be a two element list of either NoneType or float"
            )


class FluxCalibratePar(parset.ParSet):
    """
    New-style parameter set holding the arguments for how to perform the flux
    calibration (replacement for FluxCalibratePar).

    Mirrors the legacy `FluxCalibratePar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'fluxcalibrate'

    parameters = {
        'extrap_sens': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                "If False (default), the code will crash if one tries to use "
                "sensfunc at wavelengths outside its defined domain. By changing the "
                "par['sensfunc']['extrap_blu'] and par['sensfunc']['extrap_red'] this domain "
                "can be extended. If True the code will blindly extrapolate."
            ),
        ),
        'extinct_correct': parset.set_parameter_definition(
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
        'extinct_file': parset.set_parameter_definition(
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
        'use_archived_sens': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Use an archived sensfunc to flux calibration',
        ),
    }


class SensfuncUVISPar(parset.ParSet):
    """
    New-style parameter set for sensitivity function computation using the UV algorithm
    (replacement for SensfuncUVISPar).

    Mirrors the legacy `SensfuncUVISPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'sensfunc_uvis'

    parameters = {
        'std_file': parset.set_parameter_definition(
            dtype=str,
            default=None,
            descr='Standard star file to generate sensfunc',
        ),
        'std_obj_id': parset.set_parameter_definition(
            dtype=[str, int],
            default=None,
            descr=('Specifies object in spec1d file to use as standard. The brightest object found is used otherwise.'),
        ),
        'sensfunc': parset.set_parameter_definition(
            dtype=str,
            default=None,
            descr='FITS file that contains or will contain the sensitivity function.',
        ),
        'extinct_correct': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If ``extinct_correct=True`` the code will use an atmospheric extinction model to '
                'extinction correct the data below 10000A. Note that this correction makes no '
                'sense if one is telluric correcting and this shold be set to False'
            ),
        ),
        'extinct_file': parset.set_parameter_definition(
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
        'telluric_correct': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                "If ``telluric_correct=True`` the code will grab the sens_dict['telluric'] tag from the "
                "sensfunc dictionary and apply it to the data."
            ),
        ),
        'telluric': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'If ``telluric=True`` the code creates a synthetic standard star spectrum using the Kurucz models, '
                'the sens func is created setting nresln=1.5 it contains the correction for telluric lines.'
            ),
        ),
        'polycorrect': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Whether you want to correct the sensfunc with polynomial in the telluric and recombination line regions',
        ),
        'polyfunc': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Whether you want to use the polynomial fit as your final SENSFUNC',
        ),
        'nresln': parset.set_parameter_definition(
            dtype=[int, float],
            default=20,
            descr='Parameter governing the spacing of the bspline breakpoints in terms of number of resolution elements.',
        ),
        'resolution': parset.set_parameter_definition(
            dtype=[int, float],
            default=3000.0,
            descr='Expected resolution of the standard star spectrum. This should be measured from the data.',
        ),
        'trans_thresh': parset.set_parameter_definition(
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


class TelluricPar(parset.ParSet):
    """
    New-style parameter set holding telluric-correction arguments (replacement for TelluricPar).

    Mirrors the legacy `TelluricPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'telluric'

    valid_teltype = ['pca', 'grid']

    parameters = {
        'telgridfile': parset.set_parameter_definition(
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
        'tell_npca': parset.set_parameter_definition(
            dtype=int,
            default=5,
            descr='Number of telluric PCA components used. Can be set to any number from 1 to 10.',
        ),
        'teltype': parset.set_parameter_definition(
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
        'sn_clip': parset.set_parameter_definition(
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
        'resln_guess': parset.set_parameter_definition(
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
        'resln_frac_bounds': parset.set_parameter_definition(
            dtype=tuple,
            default=(0.6, 1.4),
            descr=(
                'Bounds for the resolution fit optimization which is part of the '
                'telluric model.  This range is in units of ``resln_guess``, so the '
                'default would bound the spectral resolution fit to be within the '
                'range ``bounds_resln = (0.6*resln_guess, 1.4*resln_guess)``.'
            ),
        ),
        'pix_shift_bounds': parset.set_parameter_definition(
            dtype=tuple,
            default=(-5.0, 5.0),
            descr='Bounds for the pixel shift optimization in the telluric model fit in units of pixels.  The atmosphere will be allowed to shift within this range during the fit.',
        ),
        'delta_coeff_bounds': parset.set_parameter_definition(
            dtype=tuple,
            default=(-20.0, 20.0),
            descr='Parameters setting the polynomial coefficient bounds for sensfunc optimization.',
        ),
        'minmax_coeff_bounds': parset.set_parameter_definition(
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
        'maxiter': parset.set_parameter_definition(
            dtype=int,
            default=2,
            descr='Maximum number of iterations for the telluric + object model fitting.  The code performs multiple iterations rejecting outliers at each step.  The fit is then performed anew to the remaining good pixels.  For this reason if you run with the ``disp=True`` option, you will see that the f(x) loss function gets progressively better during the iterations.',
        ),
        'sticky': parset.set_parameter_definition(
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
        'lower': parset.set_parameter_definition(
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
        'upper': parset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr='Upper rejection threshold in units of ``sigma_corr*sigma``, where ``sigma`` is the formal noise of the spectrum, and ``sigma_corr`` is an empirically determined correction to the formal error. See ``lower`` for additional detail.',
        ),
        'seed': parset.set_parameter_definition(
            dtype=int,
            default=777,
            descr='An initial seed for the differential evolution optimization, which is a random process.  The default is 777, which will be used to generate a unique seed for every order.  A specific seed is used because otherwise the random number generator will use the time for the seed, and the results will not be reproducible.',
        ),
        'tol': parset.set_parameter_definition(
            dtype=float,
            default=1e-3,
            descr='Relative tolerance for converage of the differential evolution optimization. See `scipy.optimize.differential_evolution`_ for details.',
        ),
        'popsize': parset.set_parameter_definition(
            dtype=int,
            default=30,
            descr='A multiplier for setting the total population size for the differential evolution optimization. See `scipy.optimize.differential_evolution`_ for details.',
        ),
        'recombination': parset.set_parameter_definition(
            dtype=[int, float],
            default=0.7,
            descr='The recombination constant for the differential evolution optimization. This should be in the range between 0 and 1. See `scipy.optimize.differential_evolution`_ for details.',
        ),
        'polish': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='If True then differential evolution will perform an additional optimization at the end to polish the best fit at the end, which can improve the optimization slightly. See `scipy.optimize.differential_evolution`_ for details.',
        ),
        'disp': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Argument for `scipy.optimize.differential_evolution`_ that will display status messages to the screen indicating the status of the optimization.  See documentation for :class:`~pypeit.core.telluric.Telluric` for a description of the output and how to know if things are working well.',
        ),
        'only_orders': parset.set_parameter_definition(
            dtype=[int, list, np.ndarray],
            default=None,
            descr='Order number, or list of order numbers if you only want to fit specific orders.',
        ),
        'objmodel': parset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=('The object model to be used for telluric fitting. Currently the options are: ``qso``, ``star``, and ``poly``.  For ``qso``, you might need to set ``redshift`` and ``bal_wv_min_max``.  For ``star``, you must set ``star_type``, ``star_ra``, ``star_dec``, and ``star_mag``.  For ``poly``, you might need to set ``fit_wv_min_max`` and ``norder``.'),
        ),
        'redshift': parset.set_parameter_definition(
            dtype=[int, float],
            default=0.0,
            descr='The redshift for the object model. This is currently only used by the QSO model.',
        ),
        'delta_redshift': parset.set_parameter_definition(
            dtype=float,
            default=0.1,
            descr='Range within the redshift can be varied for telluric fitting, i.e. the code performs a bounded optimization within the redshift +- delta_redshift.',
        ),
        'pca_file': parset.set_parameter_definition(
            dtype=str,
            default='qso_pca_1200_3100.fits',
            descr='Fits file containing quasar PCA model. Needed for the QSO model.  If you change the default, you might need to set ``pca_lower`` and ``pca_upper``.',
        ),
        'npca': parset.set_parameter_definition(
            dtype=int,
            default=8,
            descr='Number of pca for the objmodel=qso qso PCA fit',
        ),
        'bal_wv_min_max': parset.set_parameter_definition(
            dtype=[list, np.ndarray],
            default=None,
            descr='Min/max wavelength of broad absorption features. If there are several BAL features, the format for this mask is ``[wave_min_bal1, wave_max_bal1, wave_min_bal2, wave_max_bal2,...]``. These masked pixels will be ignored during the fitting.',
        ),
        'bounds_norm': parset.set_parameter_definition(
            dtype=tuple,
            default=(0.1, 3.0),
            descr='Normalization bounds for scaling the initial object model.',
        ),
        'tell_norm_thresh': parset.set_parameter_definition(
            dtype=[int, float],
            default=0.9,
            descr='Threshold of telluric absorption region',
        ),
        'pca_lower': parset.set_parameter_definition(
            dtype=[int, float],
            default=1220.0,
            descr='Minimum wavelength for the qso pca model',
        ),
        'pca_upper': parset.set_parameter_definition(
            dtype=[int, float],
            default=3100.0,
            descr='Maximum wavelength for the qso pca model',
        ),
        'mask_lyman_a': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Mask the blueward of Lyman-alpha line during the fitting?',
        ),
        'star_type': parset.set_parameter_definition(
            dtype=str,
            default=None,
            descr='stellar type',
        ),
        'star_mag': parset.set_parameter_definition(
            dtype=[float, int],
            default=None,
            descr='AB magnitude in V band',
        ),
        'star_ra': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr='Object right-ascension in decimal deg',
        ),
        'star_dec': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr='Object declination in decimal deg',
        ),
        'func': parset.set_parameter_definition(
            dtype=str,
            default='legendre',
            descr='Polynomial model function',
        ),
        'model': parset.set_parameter_definition(
            dtype=str,
            default='exp',
            descr='Types of polynomial model. Options are ``poly``, ``square``, ``exp`` corresponding to normal polynomial, squared polynomial, or exponentiated polynomial.',
        ),
        'polyorder': parset.set_parameter_definition(
            dtype=int,
            default=3,
            descr='Order of the polynomial model fit',
        ),
        'fit_wv_min_max': parset.set_parameter_definition(
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


class SensFuncPar(parset.ParSet):
    """
    New-style parameter set holding the arguments for sensitivity function computation
    using the UV algorithm (replacement for SensFuncPar).

    Mirrors the legacy `SensFuncPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'sensfunc'

    parameters = {
        'use_flat': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='If True, the flatfield spectrum will be used when computing the sensitivity function.',
        ),
        'extr': parset.set_parameter_definition(
            dtype=str,
            default='OPT',
            descr=(
                'Extraction method to use for the sensitivity function.  Options are: '
                "'OPT' (optimal extraction), 'BOX' (boxcar extraction). Default is 'OPT'."
            ),
        ),
        'extrap_blu': parset.set_parameter_definition(
            dtype=float,
            default=0.1,
            descr=(
                'Fraction of minimum wavelength coverage to grow the wavelength coverage of the '
                'sensitivitity function in the blue direction (`i.e.`, if the standard star spectrum '
                'cuts off at ``wave_min``) the sensfunc will be extrapolated to cover down to '
                ' (1.0 - ``extrap_blu``) * ``wave_min``'
            ),
        ),
        'extrap_red': parset.set_parameter_definition(
            dtype=float,
            default=0.1,
            descr=(
                'Fraction of maximum wavelength coverage to grow the wavelength coverage of the '
                'sensitivitity function in the red direction (`i.e.`, if the standard star spectrum'
                'cuts off at ``wave_max``) the sensfunc will be extrapolated to cover up to '
                ' (1.0 + ``extrap_red``) * ``wave_max``'
            ),
        ),
        'samp_fact': parset.set_parameter_definition(
            dtype=float,
            default=1.5,
            descr=(
                'Sampling factor to make the wavelength grid for sensitivity function finer or coarser. '
                'samp_fact > 1.0 oversamples (finer), samp_fact < 1.0 undersamples (coarser).'
            ),
        ),
        'multi_spec_det': parset.set_parameter_definition(
            dtype=list,
            default=None,
            descr=(
                'List of detectors (identified by their string name, like '
                "'DET01') to splice together for multi-detector instruments "
                '(e.g. DEIMOS). It is assumed that there is *no* overlap in '
                'wavelength across detectors (might be ok if there is).  If '
                "entered as a list of integers, they should be converted to ' "
                'the detector name.  **Cannot be used with detector mosaics.**'
            ),
        ),
        'trim_std_pixs': parset.set_parameter_definition(
            dtype=[list, tuple],
            default=None,
            descr=(
                'List or tuple of two integers specifying the number of pixels to trim'
                'from the start and end of the 1D standard star spectrum. '
                'Example: [10, 5] will trim 10 pixels from the start (blue side)'
                'and 5 pixels from the end (red side) of the spectrum. '
            ),
        ),
        'algorithm': parset.set_parameter_definition(
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
        'UVIS': parset.set_parameter_definition(
            dtype=SensfuncUVISPar,
            default=SensfuncUVISPar(),
            descr='Parameters for the UVIS sensfunc algorithm',
        ),
        'IR': parset.set_parameter_definition(
            dtype=TelluricPar,
            default=TelluricPar(),
            descr='Parameters for the IR sensfunc algorithm',
        ),
        'polyorder': parset.set_parameter_definition(
            dtype=[int, list],
            default=5,
            descr='Polynomial order for sensitivity function fitting',
        ),
        'star_type': parset.set_parameter_definition(
            dtype=str,
            default=None,
            descr='Spectral type of the standard star (for near-IR mainly)',
        ),
        'star_mag': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr='Magnitude of the standard star (for near-IR mainly)',
        ),
        'star_ra': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr='RA of the standard star. This will override values in the header (`i.e.`, if they are wrong or absent)',
        ),
        'star_dec': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr='DEC of the standard star. This will override values in the header (`i.e.`, if they are wrong or absent)',
        ),
        'mask_hydrogen_lines': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Mask hydrogen Balmer, Paschen, Brackett, and Pfund recombination lines in the sensitivity function fit. '
                'A region equal to ``hydrogen_mask_wid`` on either side of the line center is masked.'
            ),
        ),
        'hydrogen_mask_wid': parset.set_parameter_definition(
            dtype=float,
            default=10.0,
            descr='Mask width from line center for hydrogen recombination lines in Angstroms (total mask width is 2x this value).',
        ),
        'mask_helium_lines': parset.set_parameter_definition(
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
            raise ValueError(f"'extr' must be one of: {', '.join(allowed_extractions)}")

        # check trim_std_pixs format
        if self.data['trim_std_pixs'] is not None:
            if (
                not isinstance(self.data['trim_std_pixs'], (list, tuple))
                or len(self.data['trim_std_pixs']) != 2
            ):
                raise ValueError("`trim_std_pixs` must be a list or tuple of two integers.")


class SlitMaskPar(parset.ParSet):
    """
    New-style parameter set holding the arguments for slitmask ingestion and object assignment

    Mirrors the legacy `SlitMaskPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'slitmask'

    parameters = {
        'obj_toler': parset.set_parameter_definition(
            dtype=[int, float],
            default=1.0,
            descr=(
                'If slitmask design information is provided, and slit matching is performed '
                '(``use_maskdesign = True`` in ``EdgeTracePar``), this parameter provides '
                'the desired tolerance (arcsec) to match sources to targeted objects'
            ),
        ),
        'assign_obj': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='If SlitMask object was generated, assign RA,DEC,name to detected objects',
        ),
        'use_alignbox': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Use stars in alignment boxes to compute the slitmask offset. '
                'If this is set to ``True`` PypeIt will NOT compute '
                'the offset using ``snr_thrshd`` or ``bright_maskdef_id``'
            ),
        ),
        'snr_thrshd': parset.set_parameter_definition(
            dtype=[int, float],
            default=50.0,
            descr=(
                'Objects detected above this S/N threshold will '
                'be used to compute the slitmask offset. This is the default behaviour for DEIMOS '
                ' unless ``slitmask_offset``, ``bright_maskdef_id`` or ``use_alignbox`` is set.'
            ),
        ),
        'slitmask_offset': parset.set_parameter_definition(
            dtype=[int, float],
            default=None,
            descr=(
                'User-provided slitmask offset (pixels) from the position expected by '
                'the slitmask design. This is optional, and if set PypeIt will NOT compute '
                'the offset using ``snr_thrshd`` or ``bright_maskdef_id``.'
            ),
        ),
        'use_dither_offset': parset.set_parameter_definition(
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
        'bright_maskdef_id': parset.set_parameter_definition(
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
        'extract_missing_objs': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Force extraction of undetected objects at the location expected '
                  'from the slitmask design.',
        ),
        'missing_objs_fwhm': parset.set_parameter_definition(
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
        'missing_objs_boxcar_rad': parset.set_parameter_definition(
            dtype=[int, float],
            default=1.0,
            descr='Indicates the boxcar radius in arcsec for the force '
                  'extraction of undetected objects. ',
        ),
    }


class ReduxPar(parset.ParSet):
    """
    New-style parameter set for global reduction settings (replacement for ReduxPar).

    Mirrors the legacy `ReduxPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'redux'

    parameters = {
        'spectrograph': parset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=(
                'Spectrograph that provided the data to be reduced.  '
                'See :ref:`instruments` for valid options.'
            ),
        ),
        'quicklook': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Run a quick look reduction? This is usually good if you want to quickly '
                'reduce the data (usually at the telescope in real time) to get an initial '
                'estimate of the data quality.'
            ),
        ),
        'detnum': parset.set_parameter_definition(
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
        'slitspatnum': parset.set_parameter_definition(
            dtype=[str, list],
            default=None,
            descr=(
                'Restrict reduction to a set of slit DET:SPAT values (closest slit is used). '
                'Example syntax -- slitspatnum = DET01:175,DET01:205 or MSC02:2234  If you are re-running the code, '
                '(i.e. modifying one slit) you *must* have the precise SPAT_ID index.'
            ),
        ),
        'maskIDs': parset.set_parameter_definition(
            dtype=[str, int, list],
            default=None,
            descr=(
                'Restrict reduction to a set of slitmask IDs '
                'Example syntax -- ``maskIDs = 818006,818015`` '
                'This must be used with detnum (for now).'
            ),
        ),
        'sortroot': parset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=(
                'A filename given to output the details of the sorted files.  If '
                'None, the default is the root name of the pypeit file.  If off, '
                'no output is produced.'
            ),
        ),
        'calwin': parset.set_parameter_definition(
            dtype=[int, float],
            default=0,
            descr=(
                'The window of time in hours to search for calibration frames for a '
                'science frame'
            ),
        ),
        'ignore_bad_headers': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Ignore bad headers (NOT recommended unless you know it is safe).',
        ),
        'scidir': parset.set_parameter_definition(
            dtype=str,
            default='Science',
            descr='Directory relative to calling directory to write science files.',
        ),
        'qadir': parset.set_parameter_definition(
            dtype=str,
            default='QA',
            descr='Directory relative to calling directory to write quality '
                  'assessment files.',
        ),
        'redux_path': parset.set_parameter_definition(
            dtype=str,
            default=os.getcwd(),
            descr=(
                'Path to folder for performing reductions.  Default is the '
                'current working directory.'
            ),
        ),
        'chk_version': parset.set_parameter_definition(
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


class WavelengthSolutionPar(parset.ParSet):
    """
    New-style parameter set for wavelength solution settings (replacement for WavelengthSolutionPar).

    Mirrors the legacy `WavelengthSolutionPar` in :mod:`pypeit.par.pypeitpar`.
    """
    
    valid_methods = ['holy-grail', 'identify', 'reidentify', 'echelle', 'full_template']
    
    valid_reference_frames = ['observed', 'heliocentric', 'barycentric']

    parameters = {
        'reference': parset.set_parameter_definition(
            dtype=str,
            default='arc',
            options=['arc', 'sky', 'pixel'],
            descr=(
                'Perform wavelength calibration with an arc, sky frame.  Use '
                '\'pixel\' for no wavelength solution.'
            ),
        ),
        'method': parset.set_parameter_definition(
            dtype=str,
            default='holy-grail',
            options=valid_methods,
            descr=(
                'Method to use to fit the individual arc lines.  Note that some of '
                'the available methods should not be used; they are unstable and '
                'require significant parameter tweaking to succeed.  You should use '
                'one of \'holy-grail\', \'reidentify\', or \'full_template\'.  '
                '\'holy-grail\' attempts to get a first guess at line IDs by looking '
                'for patterns in the line locations.  It is fully automated.  When '
                'it works, it works well; however, it can fail catastrophically.  '
                'Instead, \'reidentify\' and \'full_template\' are the preferred '
                'methods.  They require an archived wavelength solution for your '
                'specific instrument/grating combination as a reference.  '
                'This is used to anchor the wavelength solution for the data being '
                f"reduced.  All options are: {', '.join(valid_methods)}."
            ),
        ),
        'echelle': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Is this an echelle spectrograph? If yes an additional 2-d fit '
                'wavelength fit will be performed as a function of spectral pixel '
                'and order number to improve the wavelength solution'
            ),
        ),
        'ech_2dfit': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'By default, a 2D fit to the echelle orders will be performed. If set to False, '
                'then even if this is an echelle spectrograph, the 2-d fit will not be generated. '
                'Set this to False if you wish to use the arxiv solution exactly as it '
                'was saved with pypeit_identify.'
            ),
        ),
        'ech_separate_2d': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'For echelle spectrographs, fit the 2D solutions on separate detectors separately'
            ),
        ),
        'ech_nspec_coeff': parset.set_parameter_definition(
            dtype=int,
            default=4,
            descr=(
                'For echelle spectrographs, this is the order of the final '
                '2d fit to the spectral dimension.  You should choose this '
                'to be the n_final of the fits to the individual orders.'
            ),
        ),
        'ech_norder_coeff': parset.set_parameter_definition(
            dtype=int,
            default=4,
            descr=(
                'For echelle spectrographs, this is the order of the final '
                '2d fit to the order dimension.'
            ),
        ),
        'ech_sigrej': parset.set_parameter_definition(
            dtype=[int, float],
            default=2.0,
            descr=(
                'For echelle spectrographs, this is the sigma-clipping rejection '
                'threshold in the 2d fit to spectral and order dimensions'
            ),
        ),
        'bad_orders_maxfrac': parset.set_parameter_definition(
            dtype=float,
            default=0.25,
            descr=(
                'For echelle spectrographs (i.e., ``echelle=True``), '
                'this is the maximum fraction of orders (per detector) with failed 1D fit, '
                'for PypeIt to attempt a refit.'
            ),
        ),
        'frac_rms_thresh': parset.set_parameter_definition(
            dtype=float,
            default=1.5,
            descr=(
                'For echelle spectrographs (i.e., ``echelle=True``), '
                'this is the fractional change in the RMS threshold used '
                'when a 1D fit is re-attempted for failed orders.'
            ),
        ),
        'lamps': parset.set_parameter_definition(
            dtype=list,
            default=None,
            descr=(
                'Name of one or more ions used for the wavelength calibration.  Use '
                '``None`` for no calibration. Choose ``use_header`` to use the list of lamps '
                'recorded in the header of the arc frames (this is currently '
                'available only for Keck DEIMOS, Keck LRIS, MMT Blue Channel, and LDT DeVeny).'
            ),
        ),
        'use_instr_flag': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'If True, restrict to lines matching the instrument.  WARNING: This '
                'is only implemented for shane_kast_red + HolyGrail.  Do not use it '
                'unless you really know what you are doing.'
            ),
        ),
        'sigdetect': parset.set_parameter_definition(
            dtype=[int, float, list, np.ndarray],
            default=5.0,
            descr=(
                'Sigma threshold above fluctuations for arc-line detection.  Arcs are '
                'continuum subtracted and the fluctuations are computed after continuum '
                'subtraction.  This can be a single number or a vector (list or numpy '
                'array) that provides the detection threshold for each slit.'
            ),
        ),
        'fwhm': parset.set_parameter_definition(
            dtype=[int, float],
            default=4.0,
            descr=(
                'Spectral sampling of the arc lines. This is the FWHM of an arcline in '
                'binned pixels of the input arc image. Note that this is used also in the '
                'wave tilts calibration.'
            ),
        ),
        'fwhm_fromlines': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Estimate spectral resolution in each slit using the arc lines. '
                'If True, the estimated FWHM will override ``fwhm`` in '
                'the determination of the wavelength solution (including the '
                'calculation of the threshold for the solution RMS, see '
                '``rms_thresh_frac_fwhm``), and ALSO for the wave tilts calibration.'
            ),
        ),
        'fwhm_spat_order': parset.set_parameter_definition(
            dtype=int,
            default=0,
            descr=(
                'This parameter determines the spatial polynomial order to use in the '
                '2D polynomial fit to the FWHM of the arc lines. See also, fwhm_spec_order.'
            ),
        ),
        'fwhm_spec_order': parset.set_parameter_definition(
            dtype=int,
            default=1,
            descr=(
                'This parameter determines the spectral polynomial order to use in the '
                '2D polynomial fit to the FWHM of the arc lines. See also, fwhm_spat_order.'
            ),
        ),
        'reid_arxiv': parset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=(
                'Name of the archival wavelength solution file that will be used '
                'for the wavelength reidentification.  Only used if ``method`` is '
                '\'reidentify\' or \'full_template\'.'
            ),
        ),
        'nreid_min': parset.set_parameter_definition(
            dtype=int,
            default=1,
            descr=(
                'Minimum number of times that a given candidate reidentified line must be properly '
                'matched with a line in the arxiv to be considered a good reidentification. If there '
                'is a lot of duplication in the arxiv of the spectra in question (i.e. multislit) set '
                'this to a number like 1-4. For echelle this depends on the number of solutions in the '
                'arxiv.  Set this to 1 for fixed format echelle spectrographs.  For an echelle with a '
                'tiltable grating, this will depend on the number of solutions in the arxiv.'
            ),
        ),
        'reid_cont_sub': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If True, continuum subtract the arc and arxiv spectrum before '
                'the wavelength reidentification. '
            ),
        ),
        'wvrng_arxiv': parset.set_parameter_definition(
            dtype=list,
            default=None,
            descr=(
                'Cut the arxiv template down to this specified wavelength range [min,max]'
            ),
        ),
        'nsnippet': parset.set_parameter_definition(
            dtype=int,
            default=2,
            descr=(
                "Number of spectra to chop the arc spectrum into when ``method`` is 'full_template'"
            ),
        ),
        'cc_shift_range': parset.set_parameter_definition(
            dtype=tuple,
            default=None,
            descr=(
                'Range of pixel shifts allowed when cross-correlating the '
                'input arc spectrum with the archive spectrum.  If None, '
                '``cc_offset_minmax`` will be used to determine this range.'
            ),
        ),
        'cc_thresh': parset.set_parameter_definition(
            dtype=[float, list, np.ndarray],
            default=0.70,
            descr=(
                'Threshold for the *global* cross-correlation coefficient between '
                'an input spectrum and member of the archive required to attempt '
                'reidentification.  Spectra from the archive with a lower '
                'cross-correlation are not used for reidentification. This can be '
                'a single number or a list/array providing the value for each slit.'
            ),
        ),
        'cc_local_thresh': parset.set_parameter_definition(
            dtype=float,
            default=0.70,
            descr=(
                'Threshold for the *local* cross-correlation coefficient, '
                'evaluated at each reidentified line,  between an input '
                'spectrum and the shifted and stretched archive spectrum '
                'above which a line must be to be considered a good line '
                'for reidentification. The local cross-correlation is '
                'evaluated at each candidate reidentified line (using a '
                'window of nlocal_cc), and is then used to score the the '
                'reidentified lines to arrive at the final set of good '
                'reidentifications.'
            ),
        ),
        'nlocal_cc': parset.set_parameter_definition(
            dtype=int,
            default=11,
            descr=(
                'Size of pixel window used for local cross-correlation computation for each arc line. '
                'If not an odd number one will be added to it to make it odd.'
            ),
        ),
        'rms_thresh_frac_fwhm': parset.set_parameter_definition(
            dtype=float,
            default=0.15,
            descr=(
                'Maximum RMS (expressed as fraction of the FWHM) for keeping '
                'a slit/order solution. If ``fwhm_fromlines`` is True, '
                'FWHM will be computed from the arc lines in each slits, otherwise ``fwhm`` '
                'will be used. This parameter is used for the \'holy-grail\', '
                '\'reidentify\', and \'echelle\' methods and  when re-analyzing '
                'a slit using the ``redo_slits`` parameter. '
            ),
        ),
        'match_toler': parset.set_parameter_definition(
            dtype=float,
            default=2.0,
            descr=(
                'Matching tolerance in pixels when searching for new lines. This '
                'is the difference in pixels between the wavlength assigned to '
                'an arc line by an iteration of the wavelength solution to the '
                'wavelength in the line list.  This parameter is also used as '
                'the matching tolerance in pixels for a line reidentification.  '
                'A good line match must match within this tolerance to the '
                'shifted and stretched archive spectrum, and the archive '
                'wavelength solution at this match must be within match_toler '
                'dispersion elements from the line in line list.'
            ),
        ),
        'func': parset.set_parameter_definition(
            dtype=str,
            default='legendre',
            descr=(
                'Function used for wavelength solution fits'
            ),
        ),
        'n_first': parset.set_parameter_definition(
            dtype=int,
            default=2,
            descr=(
                'Order of first guess fit to the wavelength solution.'
            ),
        ),
        'sigrej_first': parset.set_parameter_definition(
            dtype=float,
            default=2.0,
            descr=(
                'Number of sigma for rejection for the first guess to the wavelength solution.'
            ),
        ),
        'n_final': parset.set_parameter_definition(
            dtype=[int, float, list, np.ndarray],
            default=4,
            descr=(
                'Order of final fit to the wavelength solution (there are n_final+1 parameters '
                'in the fit). This can be a single number or a list/array providing the value '
                'for each slit'
            ),
        ),
        'sigrej_final': parset.set_parameter_definition(
            dtype=float,
            default=3.0,
            descr=(
                'Number of sigma for rejection for the final guess to the wavelength solution.'
            ),
        ),
        'numsearch': parset.set_parameter_definition(
            dtype=int,
            default=20,
            descr=(
                'Number of brightest arc lines to search for in preliminary identification'
            ),
        ),
        'nfitpix': parset.set_parameter_definition(
            dtype=int,
            default=5,
            descr=(
                'Number of pixels to fit when deriving the centroid of the arc lines (an odd number is best)'
            ),
        ),
        'refframe': parset.set_parameter_definition(
            dtype=str,
            default='heliocentric',
            options=valid_reference_frames,
            descr=(
                'Frame of reference for the wavelength calibration.  Options are: '
                f'{", ".join(valid_reference_frames)}'
            )
        ),
        'redo_slits': parset.set_parameter_definition(
            dtype=[int, list],
            descr=(
                'Redo the input slit(s) [multislit] or order(s) [echelle]'
            ),
        ),
        'qa_log': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Governs whether the wavelength solution arc line QA plots will have log or '
                'linear scaling.  If True, the scaling will be log, if False linear'
            ),
        ),
        'cc_percent_ceil': parset.set_parameter_definition(
            dtype=float,
            default=50.0,
            descr=(
                'Determines the percentile at which to cap lines used in cross correlation, '
                'to prevent large lines from dominating. If 100, all lines are allowed at their '
                'maximum heights. May produce spurious peaks in xcorr'
            ),
        ),
        'echelle_pad': parset.set_parameter_definition(
            dtype=int,
            default=3,
            descr=(
                'Number of orders by which to pad the echellogram reference in the echelle '
                'method. Values > 0 allow for some error in the reddest order guess, but '
                'require sufficient reference orders.'
            ),
        ),
        'cc_offset_minmax': parset.set_parameter_definition(
            dtype=float,
            default=1.0,
            descr=(
                'Fraction of the total spectral pixels used to determine the range of '
                'pixel shifts allowed when cross-correlating the input arc spectrum with '
                'the archive spectrum. Restricting this can be crucial if there are few '
                'reference lines and the cross correlation can get confused. '
                'This parameter is only used if ``cc_shift_range`` is None.'
            ),
        ),
        'stretch_func': parset.set_parameter_definition(
            dtype=str,
            default='quadratic',
            options=['linear', 'quadratic'],
            descr=(
                'Whether to use a linear (linear) or quadratic (quad) function to stretch '
                'the extracted arcs when identifying emission lines with reidentify. For '
                'NIRSPEC, the quadratic mode tends to do better because the wavelength '
                'solution is typically at least 2nd or 3rd order.'
            ),
        ),
    }

    


class EdgeTracePar(parset.ParSet):
    """
    New-style parameter set for slit edge tracing (replacement for EdgeTracePar).

    Mirrors the legacy `EdgeTracePar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'edgetrace'

    card_prefix = 'ETP'

    parameters = {
        'filt_iter': parset.set_parameter_definition(
            dtype=int,
            default=0,
            descr=(
                'Number of median-filtering iterations to perform on sqrt(trace) '
                'image before applying to Sobel filter to detect slit/order edges.'
            ),
        ),
        'sobel_mode': parset.set_parameter_definition(
            dtype=str,
            default='nearest',
            options=['nearest', 'constant'],
            descr=(
                'Mode for Sobel filtering.  Default is \'nearest\'; note we find'
                '\'constant\' works best for DEIMOS.'
            ),
        ),
        'edge_thresh': parset.set_parameter_definition(
            dtype=[int, float],
            default=20.0,
            descr=(
                'Threshold for finding edges in the Sobel-filtered significance image.'
            ),
        ),
        'sobel_enhance': parset.set_parameter_definition(
            dtype=int,
            default=0,
            descr=(
                'Enhance the sobel filtering? A value of 0 will not enhance the sobel filtering. '
                'Any other value > 0 will sum the sobel values. For example, a value of 3 will '
                'combine the sobel values for the 3 nearest pixels. This is useful when a slit '
                'edge is poorly defined (e.g. vignetted).'
            ),
        ),
        'exclude_regions': parset.set_parameter_definition(
            dtype=[list, str],
            default=None,
            descr=(
                'User-defined regions to exclude from the slit tracing. To set this parameter, '
                'the text should be a comma separated list of pixel ranges (in the x direction) '
                'to be excluded and the detector number. For example, the following string '
                '1:0:20,1:300:400  would select two regions in det=1 between pixels 0 and 20 '
                'and between 300 and 400.'
            ),
        ),
        'follow_span': parset.set_parameter_definition(
            dtype=int,
            default=20,
            descr=(
                'In the initial connection of spectrally adjacent edge '
                'detections, this sets the number of previous spectral rows '
                'to consider when following slits forward.'
            ),
        ),
        'det_min_spec_length': parset.set_parameter_definition(
            dtype=[int, float],
            default=0.33,
            descr=(
                'The minimum spectral length (as a fraction of the ' 
                'detector size) of a trace determined by direct ' 
                'measurements of the detector data (as opposed to what ' 
                'should be included in any modeling approach; see ' 
                'fit_min_spec_length).'
            ),
        ),
        'trim_spec': parset.set_parameter_definition(
            dtype=list,
            default=None,
            descr=(
                'User-defined truncation of all slits in the spectral direction.'
                'Should be two integers, e.g. 100,150 trims 100 pixels from the '
                'short wavelength end and 150 pixels from the long wavelength '
                'end of the spectral axis of the detector.'
            ),
        ),
        'mask_off_detector': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Mask spectral regions in each slit/order where more than '
                '50% of the slit spatial coverage falls off the detector. '
            ),
        ),
        'max_shift_abs': parset.set_parameter_definition(
            dtype=[int, float],
            default=0.5,
            descr=(
                'Maximum spatial shift in pixels between an input edge '
                'location and the recentroided value.'
            ),
        ),
        'max_shift_adj': parset.set_parameter_definition(
            dtype=[int, float],
            default=0.15,
            descr=(
                'Maximum spatial shift in pixels between the edges in '
                'adjacent spectral positions.'
            ),
        ),
        'max_spat_error': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Maximum error in the spatial position of edges in pixels.'
            ),
        ),
        'match_tol': parset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr=(
                'Same-side slit edges below this separation in pixels are '
                'considered part of the same edge.'
            ),
        ),
        'fit_function': parset.set_parameter_definition(
            dtype=str,
            default='legendre',
            options=['polynomial', 'legendre', 'chebyshev'],
            descr=(
                'Function fit to edge measurements.  Options are: polynomial, legendre, chebyshev'
            ),
        ),
        'fit_order': parset.set_parameter_definition(
            dtype=int,
            default=5,
            descr=(
                'Order of the function fit to edge measurements.'
            ),
        ),
        'fit_maxdev': parset.set_parameter_definition(
            dtype=[int, float],
            default=5.0,
            descr=(
                'Maximum deviation between the fitted and measured edge position '
                'for rejection in spatial pixels.'
            ),
        ),
        'fit_maxiter': parset.set_parameter_definition(
            dtype=int,
            default=25,
            descr=(
                'Maximum number of rejection iterations during edge fitting.'
            ),
        ),
        'fit_niter': parset.set_parameter_definition(
            dtype=int,
            default=1,
            descr=(
                'Number of iterations of re-measuring and re-fitting the edge '
                'data; see :func:`~pypeit.core.trace.fit_trace`.'
            ),
        ),
        'fit_min_spec_length': parset.set_parameter_definition(
            dtype=float,
            default=0.6,
            descr=(
                'Minimum unmasked spectral length of a traced slit edge '
                'to use in any modeling procedure (polynomial fitting '
                'or PCA decomposition).'
            ),
        ),
        'auto_pca': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'During automated tracing, attempt to construct a PCA decomposition '
                'of the traces. When True, the edge traces resulting from the '
                'initial detection, centroid refinement, and polynomial fitting '
                'must meet a set of criteria for performing the pca; see ' 
                ':func:`pypeit.edgetrace.EdgeTraceSet.can_pca`.  If False, the '
                '``sync_predict`` parameter *cannot* be set to ``pca``; if it is '
                'not, the value is set to ``nearest`` and a warning is issued when '
                'validating the parameter set.'
            ),
        ),
        'left_right_pca': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Construct a PCA decomposition for the left and right traces '
                'separately.  This can be important for cross-dispersed '
                'echelle spectrographs (e.g., Keck-NIRES)'
            ),
        ),
        'pca_min_edges': parset.set_parameter_definition(
            dtype=int,
            default=4,
            descr=(
                'Minimum number of edge traces required to perform a PCA '
                'decomposition of the trace form.  If left_right_pca is True, '
                'this minimum applies to the number of left and right traces '
                'separately.'
            ),
        ),
        'pca_n': parset.set_parameter_definition(
            dtype=int,
            descr=(
                'The number of PCA components to keep, which must be less than the '
                'number of detected traces.  If not provided, determined by '
                'calculating the minimum number of components required to explain a '
                'given percentage of variance in the edge data; see `pca_var_percent`.'
            ),
        ),
        'pca_var_percent': parset.set_parameter_definition(
            dtype=[int, float],
            default=99.8,
            descr=(
                'The percentage (i.e., not the fraction) of the variance in '
                'the edge data accounted for by the PCA used to truncate '
                'the number of PCA coefficients to keep (see `pca_n`).  '
                'Ignored if `pca_n` is provided directly.'
            ),
        ),
        'pca_function': parset.set_parameter_definition(
            dtype=str,
            default='polynomial',
            options=['polynomial', 'legendre', 'chebyshev'],
            descr=(
                'Type of function fit to the PCA coefficients for each '
                'component.  Options are: polynomial, legendre, chebyshev'
            ),
        ),
        'pca_order': parset.set_parameter_definition(
            dtype=int,
            default=2,
            descr=(
                'Order of the function fit to the PCA coefficients.'
            ),
        ),
        'pca_sigrej': parset.set_parameter_definition(
            dtype=[int, float, list],
            default=[2.0, 2.0],
            descr=(
                'Sigma rejection threshold for fitting PCA components. Individual '
                'numbers are used for both lower and upper rejection. A list of '
                'two numbers sets these explicitly (e.g., [2., 3.]).'
            ),
        ),
        'pca_maxrej': parset.set_parameter_definition(
            dtype=int,
            default=1,
            descr=(
                'Maximum number of PCA coefficients rejected during a given fit '
                'iteration.'
            ),
        ),
        'pca_maxiter': parset.set_parameter_definition(
            dtype=int,
            default=25,
            descr=(
                'Maximum number of rejection iterations when fitting the PCA '
                'coefficients.'
            ),
        ),
        'smash_range': parset.set_parameter_definition(
            dtype=list,
            default=[0.0, 1.0],
            descr=(
                'Range of the slit in the spectral direction (in fractional '
                'units) to smash when searching for slit edges.  If the '
                'spectrum covers only a portion of the image, use that range.'
            ),
        ),
        'edge_detect_clip': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Sigma clipping level for peaks detected in the collapsed, '
                'Sobel-filtered significance image.'
            ),
        ),
        'trace_median_frac': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'After detection of peaks in the rectified Sobel-filtered ' 
                'image and before refitting the edge traces, the rectified ' 
                'image is median filtered with a kernel width of ' 
                '`trace_median_frac*nspec` along the spectral dimension.'
            ),
        ),
        'trace_thresh': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'After rectification and median filtering of the Sobel-filtered ' 
                'image (see `trace_median_frac`), values in the median-filtered ' 
                'image *below* this threshold are masked in the refitting of ' 
                'the edge trace data.  If None, no masking applied.'
            ),
        ),
        'trace_rms_tol': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'After retracing edges using peaks detected in the rectified ' 
                'and collapsed image, the RMS difference (in pixels) between ' 
                'the original and refit traces are calculated.  This sets the ' 
                'upper limit of the RMS for traces that will be removed.  If ' 
                'None, no limit is set and all new traces are kept.'
            ),
        ),
        'fwhm_uniform': parset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr=(
                'The `fwhm` parameter to use when using uniform weighting in ' 
                ':func:`~pypeit.core.trace.fit_trace` when refining the PCA ' 
                'predictions of edges.  See description of ' 
                ':func:`~pypeit.core.trace.peak_trace`.'
            ),
        ),
        'niter_uniform': parset.set_parameter_definition(
            dtype=int,
            default=9,
            descr=(
                'The number of iterations of ' 
                ':func:`~pypeit.core.trace.fit_trace` to use when using ' 
                'uniform weighting.'
            ),
        ),
        'fwhm_gaussian': parset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr=(
                'The `fwhm` parameter to use when using Gaussian weighting in ' 
                ':func:`~pypeit.core.trace.fit_trace` when refining the PCA ' 
                'predictions of edges.  See description of ' 
                ':func:`~pypeit.core.trace.peak_trace`.'
            ),
        ),
        'niter_gaussian': parset.set_parameter_definition(
            dtype=int,
            default=6,
            descr=(
                'The number of iterations of ' 
                ':func:`~pypeit.core.trace.fit_trace` to use when using ' 
                'Gaussian weighting.'
            ),
        ),
        'min_edge_side_sep': parset.set_parameter_definition(
            dtype=[int, float],
            default=5.0,
            descr=(
                'Minimum separation between same-side edges (e.g., the ' 
                'minimum separation between two subsequent right-edge ' 
                'detections) in units of ``fwhm_gaussian``.  For example, ' 
                'if ``fwhm_gaussian = 3.0`` and ``min_edge_sid_sep = 5.``, ' 
                'the separation between subsequent right edges must be at ' 
                'least 15 pixels.'
            ),
        ),
        'det_buffer': parset.set_parameter_definition(
            dtype=int,
            default=5,
            descr=(
                'The minimum separation between the detector edges and a slit ' 
                'edge for any added edge traces.  Must be positive.'
            ),
        ),
        'max_nudge': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'If parts of any (predicted) trace fall off the detector edge, '
                'allow them to be nudged away from the detector edge up to and '
                'including this maximum number of pixels.  If None, no limit is '
                'set; otherwise should be 0 or larger.'
            ),
        ),
        'sync_predict': parset.set_parameter_definition(
            dtype=str,
            default='pca',
            options=['pca', 'nearest', 'auto'],
            descr=(
                'Mode to use when predicting the form of the trace to insert.  '
                'Use `pca` to use the PCA decomposition, `nearest` to '
                'reproduce the shape of the nearest trace, or `auto` to let PypeIt '
                'decide which mode to use between `pca` and `nearest`. In general, '
                'it will first try `pca`, and if that is not possible, it will use `nearest`.'
            ),
        ),
        'sync_center': parset.set_parameter_definition(
            dtype=str,
            default='median',
            options=['median', 'nearest', 'gap'],
            descr=(
                'Mode to use for determining the location of traces to insert.  '
                'Use `median` to use the median of the matched left and right '
                'edge pairs, `nearest` to use the length of the nearest slit, '
                'or `gap` to offset by a fixed gap width from the next slit edge.'
            ),
        ),
        'gap_offset': parset.set_parameter_definition(
            dtype=[int, float],
            default=5.0,
            descr=(
                'Offset (pixels) used for the slit edge gap width when inserting '
                'slit edges (see `sync_center`) or when nudging predicted slit '
                'edges to avoid slit overlaps.  This should be larger than '
                '`minimum_slit_gap` when converted to arcseconds.'
            ),
        ),
        'sync_to_edge': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If adding a first left edge or a last right edge, ignore '
                '`center_mode` for these edges and place them at the edge of '
                'the detector (with the relevant shape).'
            ),
        ),
        'bound_detector': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'When the code is ready to synchronize the left/right trace '
                'edges, the traces should have been constructed, vetted, and '
                'cleaned. This can sometimes lead to *no* valid traces. This '
                'parameter dictates what to do next. If ``bound_detector`` is '
                'True, the code will artificially add left and right edges '
                'that bound the detector; if False, the code identifies the '
                'slit-edge tracing as being unsuccessful, warns the user, and '
                'ends gracefully. Note that setting ``bound_detector`` to '
                'True is needed for some long-slit data where the slit '
                'edges are, in fact, beyond the edges of the detector.'
            ),
        ),
        'minimum_slit_dlength': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Minimum *change* in the slit length (arcsec) as a '
                'function of wavelength in arcsec.  This is mostly '
                'meant to catch cases when the polynomial fit to the '
                'detected edges becomes ill-conditioned (e.g., when '
                'the slits run off the edge of the detector) and leads '
                'to wild traces.  If reducing the order of the '
                'polynomial (``fit_order``) does not help, try using '
                'this to remove poorly constrained slits.'
            ),
        ),
        'dlength_range': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Similar to ``minimum_slit_dlength``, but constrains the '
                '*fractional* change in the slit length as a function of '
                'wavelength.  For example, a value of 0.2 means that slit '
                'length should not vary more than 20%'
                'as a function of wavelength.  '
            ),
        ),
        'minimum_slit_length': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Minimum slit length in arcsec.  Slit lengths are '
                'determined by the median difference between the left '
                'and right edge locations for the unmasked trace '
                'locations.  This is used to identify traces that are '
                '*erroneously* matched together to form slits.  Short '
                'slits are expected to be ignored or removed (see '
                ' ``clip``).  If None, no minimum slit length applied.'
            ),
        ),
        'minimum_slit_length_sci': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Minimum slit length in arcsec for a science slit.  '
                'Slit lengths are determined by the median difference '
                'between the left and right edge locations for the '
                'unmasked trace locations.  Used in combination with '
                '``minimum_slit_length``, this parameter is used to '
                'identify box or alignment slits; i.e., those slits '
                'that are shorter than ``minimum_slit_length_sci`` but '
                'larger than ``minimum_slit_length`` are box/alignment '
                'slits.  Box slits are *never* removed (see ``clip``), '
                'but no spectra are extracted from them.  If None, no '
                'minimum science slit length is applied.'
            ),
        ),
        'length_range': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Allowed range in slit length compared to the median slit '
                'length.  For example, a value of 0.3 means that slit lengths '
                'should not vary more than 30%.  Relatively shorter or longer '
                'slits are masked or clipped.  Most useful for echelle or '
                'multi-slit data where the slits should have similar or '
                'identical lengths.'
            ),
        ),
        'minimum_slit_gap': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Minimum slit gap in arcsec.  Gaps between slits are '
                'determined by the median difference between the right '
                'and left edge locations of adjacent slits.  Slits with '
                'small gaps are merged by removing the intervening traces.'
                'If None, no minimum slit gap is applied.  This should be '
                'smaller than `gap_offset` when converted to pixels.'
            ),
        ),
        'clip': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Remove traces flagged as bad, instead of only masking them.  This '
                'is currently only used by ' 
                ':func:`~pypeit.edgetrace.EdgeTraceSet.centroid_refine`.'
            ),
        ),
        'order_match': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Orders for *fixed-format* echelle spectrographs are always '
                'matched to a predefined expectation for the number of orders '
                'found and their relative placement in the detector.  This sets '
                'the tolerance allowed for matching identified "slits" to '
                'echelle orders. Must be relative to the fraction of the '
                'detector spatial scale (i.e., a value of 0.05 means that the '
                'order locations must be within 5% of the expected value).  If '
                'None, no limit is used.'
            ),
        ),
        'order_offset': parset.set_parameter_definition(
            dtype=[int, float],
            descr=(
                'Orders for *fixed-format* echelle spectrographs are always '
                'matched to a predefined expectation for the number of orders '
                'found and their relative placement in the detector.  This sets '
                'the offset to introduce to the expected order positions to '
                'improve the match for this specific data. This is an additive '
                'offset to the measured slit positions; i.e., this should '
                'minimize the difference between the expected order positions '
                'and ``self.slit_spatial_center() + offset``. Must be in the '
                'fraction of the detector spatial scale. If None, no offset '
                'is applied.'
            ),
        ),
        'add_missed_orders': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'For any Echelle spectrograph (fixed-format or otherwise), '
                'attempt to add orders that have been missed by the '
                'automated edge tracing algorithm.  For *fixed-format* '
                'echelles, this is based on the expected positions on '
                'on the detector.  Otherwise, the detected orders are '
                'modeled and used to predict the locations of missed '
                'orders; see additional parameters '
                '``order_width_poly``, ``order_gap_poly``, '
                '``order_fitrej``, ``order_outlier``, and '
                '``order_spat_range``.'
            ),
        ),
        'order_width_poly': parset.set_parameter_definition(
            dtype=int,
            default=2,
            descr=(
                'Order of the Legendre polynomial used to model the '
                'spatial width of each order as a function of spatial '
                'pixel position.  See ``add_missed_orders``.'
            ),
        ),
        'order_gap_poly': parset.set_parameter_definition(
            dtype=int,
            default=3,
            descr=(
                'Order of the Legendre polynomial used to model the spatial '
                'gap between orders as a function of the order spatial '
                'position.  See ``add_missed_orders``.'
            ),
        ),
        'order_fitrej': parset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr=(
                'When fitting the width of and gap beteween echelle orders with '
                'Legendre polynomials, this is the sigma-clipping threshold '
                'when excluding data from the fit.  See ``add_missed_orders``.'
            ),
        ),
        'order_outlier': parset.set_parameter_definition(
            dtype=[int, float],
            default=None,
            descr=(
                'When fitting the width of echelle orders with Legendre '
                'polynomials, this is the sigma-clipping threshold used to '
                'identify outliers.  Orders clipped by this threshold are '
                '*removed* from further consideration, whereas orders clipped '
                'by ``order_fitrej`` are excluded from the polynomial fit '
                'but are not removed.  Note this is *only applied to the order '
                'widths*, not the order gaps.  If None, no "outliers" are '
                'identified/removed.  Should be larger or equal to '
                '``order_fitrej``.'
            ),
        ),
        'order_spat_range': parset.set_parameter_definition(
            dtype=list,
            descr=(
                'The spatial range of the detector/mosaic over which to '
                'predict order locations.  If None, the full '
                'detector/mosaic range is used.  See ``add_missed_orders``.'
            ),
        ),
        'overlap': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Assume slits identified as abnormally short are actually due to '
                'overlaps between adjacent slits/orders.  If set to True, you *must* '
                'have also used ``length_range`` to identify left-right edge pairs '
                'that have an abnormally short separation.  For those short slits, '
                'the code attempts to convert the short slits into slit gaps.  This '
                'is particularly useful for blue orders in Keck-HIRES data.'
            ),
        ),
        'max_overlap': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                'When adding missing echelle orders based on where existing '
                'orders are found, the prediction can yield overlapping orders.  '
                'The edges of these orders are adjusted to eliminate the '
                'overlap, and orders can be added up over the spatial range of '
                'the detector set by ``order_spate_range``.  If this value is '
                'None, orders are added regardless of how much they overlap.  '
                'If not None, this defines the maximum fraction of an order '
                'spatial width that can overlap with other orders.  For example, '
                'if ``max_overlap=0.5``, any order that overlaps its neighboring '
                'orders by more than 50% will not be added as a missing order.'
            ),
        ),
        'use_maskdesign': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Use slit-mask designs to identify slits.'
            ),
        ),
        'maskdesign_filename': parset.set_parameter_definition(
            dtype=[str, list],
            default=None,
            descr=(
                'Mask design info contained in this file or files (comma separated)'
            ),
        ),
        'maskdesign_maxsep': parset.set_parameter_definition(
            dtype=[int, float],
            default=50,
            descr=(
                'Maximum allowed offset in pixels between the slit edges '
                'defined by the slit-mask design and the traced edges.'
            ),
        ),
        'maskdesign_step': parset.set_parameter_definition(
            dtype=[int, float],
            default=1,
            descr=(
                'Step in pixels used to generate a list of possible offsets '
                '(within +/- `maskdesign_maxsep`) between the slit edges defined '
                'by the mask design and the traced edges.'
            ),
        ),
        'maskdesign_sigrej': parset.set_parameter_definition(
            dtype=[int, float],
            default=3,
            descr=(
                'Number of sigma for sigma-clipping rejection during slit-mask '
                'design matching.'
            ),
        ),
        'maskdesign_trim': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'If True, the mask design information is used to trim each '
                'slit in the spectral direction. This functionality is '
                'only used for spectrographs with slit-mask designs that '
                'have information on the spectral extent of each slit (currently, '
                'only Gemini GMOS N/S).'
            ),
        ),
        'maskdesign_trim_shift': parset.set_parameter_definition(
            dtype=[int, float],
            default=0,
            descr=(
                'Shift in pixels to apply to the mask design information '
                'when trimming the slits in the spectral direction.  This '
                'is useful for cases where the mask design information '
                'is not perfectly aligned with the detector.  This '
                'functionality is only used for spectrographs with '
                'slit-mask designs that have information on the spectral '
                'extent of each slit (currently, only Gemini GMOS N/S).'
            ),
        ),
        'pad': parset.set_parameter_definition(
            dtype=int,
            default=0,
            descr=(
                'Integer number of pixels to consider beyond the slit edges when '
                "selecting pixels that are 'on' the slit."
            ),
        ),
        'add_slits': parset.set_parameter_definition(
            dtype=[str, list],
            descr=(
                'Add one or more user-defined slits.  The syntax to define a '
                'slit to add is: \'det:spec:spat_left:spat_right\' where '
                'det=detector, spec=spectral pixel, spat_left=spatial pixel of '
                'left slit boundary, and spat_righ=spatial pixel of right slit '
                'boundary.  **Multiple entries must be separated by a semi-colon.** '
                'For example, \'2:2000:2121:2322; 3:2000:1201:1500\' will add a '
                'slit to detector 2 passing through spec=2000 extending spatially '
                'from 2121 to 2322 and another on detector 3 at spec=2000 '
                'extending from 1201 to 1500.  For mosaics, use the tuple '
                'definition of the mosaic.  For example, '
                '\'(1,2,3):1537:297.2:353.5\', adds a slit that passes through '
                '(1537,297.2) on the left and (1537,353.5) on the right in the '
                'mosaic made up of detectors 1, 2, and 3.'
            ),
        ),
        'add_predict': parset.set_parameter_definition(
            dtype=str,
            default='nearest',
            descr=(
                'Sets the method used to predict the shape of the left and right ' 
                'traces for a user-defined slit inserted.  Options are (1) ' 
                '``straight`` inserts traces with a constant spatial pixels ' 
                'position, (2) ``nearest`` inserts traces with a form identical ' 
                'to the automatically identified trace at the nearest spatial ' 
                'position to the inserted slit, or (3) ``pca`` uses the PCA ' 
                'decomposition to predict the shape of the traces.'
            ),
        ),
        'rm_slits': parset.set_parameter_definition(
            dtype=[str, list],
            descr=(
                'Remove one or more user-specified slits.  The syntax used to ' 
                'define a slit to remove is: \'det:spec:spat\' where det=detector, ' 
                'spec=spectral pixel, spat=spatial pixel.  **Multiple entries must ' 
                'be separated by a semi-colon.**  For example, ' 
                "'2:2000:2121; 3:2000:1500' will remove the slit on detector 2 " 
                'that contains pixel (spec,spat)=(2000,2121) and on detector 3 ' 
                'that contains pixel (2000,1500).  For mosaics, use the tuple ' 
                'definition of the mosaic.  For example \'(1,2,3):1500:331\', ' 
                'removes the slit that contains pixel (1500,331) in the mosaic made ' 
                'up of detectors 1, 2, and 3.'
            ),
        ),
    }

    def validate(self):
        # Check if the user defined a slit to remove in a mosaic
        for k in ['rm_slits', 'add_slits']:
            if self.data.get(k) is None:
                continue
            self.data[k] = ';'.join(parse.fix_config_par_image_location(self.data[k]))

        if not self.data.get('auto_pca', True) and self.data.get('sync_predict') == 'pca':
            import warnings
            warnings.warn('sync_predict cannot be pca if auto_pca is False.  Setting to nearest.')
            self.data['sync_predict'] = 'nearest'

        if self.data.get('max_overlap') is not None and (self.data['max_overlap'] < 0 or self.data['max_overlap'] > 1):
            raise ValueError('If defined, max_overlap must be in the range [0,1].')

        if self.data.get('order_outlier') is not None and self.data.get('order_outlier') < self.data.get('order_fitrej', 0):
            log.warning('Order outlier threshold should not be less than the rejection threshold.')



class WaveTiltsPar(parset.ParSet):
    """
    New-style parameter set for tracing the monochromatic tilt along the slit

    Mirrors the legacy `WaveTiltsPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'wavetilts'

    parameters = {
        'idsonly': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Only use the arc lines that have an identified wavelength to trace '
                'tilts (CURRENTLY NOT USED!)'
            ),
        ),
        'tracethresh': parset.set_parameter_definition(
            dtype=[int, float, list, np.ndarray],
            default=20.0,
            descr=(
                'Significance threshold for arcs to be used in tracing wavelength tilts. '
                'This can be a single number or a list/array providing the value for each slit/order.'
            ),
        ),
        'sig_neigh': parset.set_parameter_definition(
            dtype=[int, float],
            default=10.0,
            descr=(
                'Significance threshold for arcs to be used in line identification for the purpose of identifying neighboring lines. '
                'The tracethresh parameter above determines the significance threshold of lines that will be traced, but these lines '
                ' must be at least nfwhm_neigh fwhm away from neighboring lines. This parameter determines the significance above which '
                ' a line must be to be considered a possible colliding neighbor. A low value of sig_neigh will result in an overall '
                ' larger number of lines, which will result in more lines above tracethresh getting rejected'
            ),
        ),
        'nfwhm_neigh': parset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr=(
                'Required separation between neighboring arc lines for them to be considered for tilt tracing in units of the '
                'the spectral fwhm (see wavelength parset where fwhm is defined)'
            ),
        ),
        'maxdev_tracefit': parset.set_parameter_definition(
            dtype=[int, float],
            default=0.2,
            descr=(
                'Maximum absolute deviation (in units of fwhm) for the legendre polynomial fits to individual '
                'arc line tilt fits during iterative trace fitting (flux weighted, then gaussian weighted)'
            ),
        ),
        'sigrej_trace': parset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr=(
                'Outlier rejection significance to determine which traced arc lines should be included in the global fit'
            ),
        ),
        'spat_order': parset.set_parameter_definition(
            dtype=[int, float, list, np.ndarray],
            default=3,
            descr=(
                'Order of the legendre polynomial to be fit to the tilt of an arc line. This parameter determines '
                'both the order of the *individual* arc line tilts, as well as the order of the spatial direction of the '
                '2d legendre polynomial (spatial, spectral) that is fit to obtain a global solution for the tilts across the '
                'slit/order. This can be a single number or a list/array providing the value for each slit'
            ),
        ),
        'spec_order': parset.set_parameter_definition(
            dtype=[int, float, list, np.ndarray],
            default=4,
            descr=(
                'Order of the spectral direction of the 2d legendre polynomial (spatial, spectral) that is '
                'fit to obtain a global solution for the tilts across the slit/order. '
                'This can be a single number or a list/array providing the value for each slit'
            ),
        ),
        'minmax_extrap': parset.set_parameter_definition(
            dtype=[list, np.ndarray],
            default=[150.0, 1000.0],
            descr=(
                'Sets how far below the last measured tilt line is extrapolated in tracewave.fit_tilts()'
            ),
        ),
        'func2d': parset.set_parameter_definition(
            dtype=str,
            default='legendre2d',
            descr=(
                'Type of function for 2D fit'
            ),
        ),
        'maxdev2d': parset.set_parameter_definition(
            dtype=[int, float],
            default=0.25,
            descr=(
                'Maximum absolute deviation (in units of fwhm) rejection threshold used to determines which pixels in global 2d fits to '
                'arc line tilts are rejected because they deviate from the model by more than this value'
            ),
        ),
        'sigrej2d': parset.set_parameter_definition(
            dtype=[int, float],
            default=3.0,
            descr=(
                'Outlier rejection significance determining which pixels on a fit to an arc line tilt '
                'are rejected by the global 2D fit'
            ),
        ),
        'rm_continuum': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Before tracing the line center at each spatial position, '
                'remove any low-order continuum in the 2D spectra.'
            ),
        ),
        'cont_rej': parset.set_parameter_definition(
            dtype=[int, float, list, np.ndarray],
            default=[3, 1.5],
            descr=(
                'The sigma threshold for rejection.  Can be a single number or two '
                'numbers that give the low and high sigma rejection, respectively.'
            ),
        ),
    }

    def validate(self):
        if hasattr(self.data['cont_rej'], '__len__'):
            if len(self.data['cont_rej']) != 2:
                raise ValueError('Continuum rejection threshold must be a single number or a '
                                 'two-element list/array.')


class FindObjPar(parset.ParSet):
    """
    New-style parameter set for finding and tracing objects (replacement for FindObjPar).

    Mirrors the legacy `FindObjPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'findobj'

    parameters = {
        'trace_npoly': parset.set_parameter_definition(
            dtype=int,
            default=5,
            descr='Order of legendre polynomial fits to object traces.',
        ),
        'maxnumber_sci': parset.set_parameter_definition(
            dtype=int,
            default=10,
            descr=(
                'Maximum number of objects to extract in a science frame.  Use '
                'None for no limit. This parameter can be useful in situations where systematics lead to '
                'spurious extra objects. Setting this parameter means they will be trimmed. '
                'For mulitslit maxnumber applies per slit, for echelle observations this '
                'applies per order. Note that objects on a slit/order impact the sky-modeling and so '
                'maxnumber should never be lower than the true number of detectable objects on your slit. '
                'For image differenced observations with positive and negative object traces, maxnumber applies '
                'to the number of positive (or negative) traces individually. In other words, if you had two positive objects and '
                'one negative object, then you would set maxnumber to be equal to two (not three). Note that if manually '
                'extracted apertures are explicitly requested, they do not count against this maxnumber. If more than '
                'maxnumber objects are detected, then highest S/N ratio objects will be the ones that are kept. '
                'For multislit observations the choice here depends on the slit length. For echelle observations '
                'with short slits we set the default to be 1'
            ),
        ),
        'maxnumber_std': parset.set_parameter_definition(
            dtype=int,
            default=5,
            descr=(
                'Maximum number of objects to extract in a standard star frame.  Same functionality as '
                'maxnumber_sci documented above. For multislit observations the default here is 5, for echelle '
                'observations the default is 1'
            ),
        ),
        'snr_thresh': parset.set_parameter_definition(
            dtype=[int, float],
            default=10.0,
            descr='S/N threshold for object finding in wavelength direction smashed image.',
        ),
        'find_trim_edge': parset.set_parameter_definition(
            dtype=list,
            default=[5, 5],
            descr='Trim the slit by this number of pixels left/right before finding objects',
        ),
        'trace_extrap_npoly': parset.set_parameter_definition(
            dtype=int,
            default=3,
            descr=(
                'Polynomial order used for trace extrapolation.  NOTE: Not consumed by the code at present. (For ``pypeit<=1.18.x``, this '
                'parameter was called ``find_extrap_npoly``.)'
            ),
        ),
        'trace_maxdev': parset.set_parameter_definition(
            dtype=[int, float],
            default=2.0,
            descr=(
                'Maximum deviation of pixels from polynomial fit to trace used to reject bad pixels in trace fitting.  (For ``pypeit<=1.18.x``, this '
                'parameter was called ``find_maxdev``.)'
            ),
        ),
        'trace_maxshift': parset.set_parameter_definition(
            dtype=[int, float],
            default=1.0,
            descr=(
                'Maximum shift allowed between the input and recalculated centroid in trace fitting.  This parameter may be increased to '
                'allow the fiter to follow curved traces (*e.g.*, for wide spectral ranges at high airmass).'
            ),
        ),
        'trace_min_max': parset.set_parameter_definition(
            dtype=list,
            default=None,
            descr=(
                'It defines the minimum and maximum pixel in the spectral direction with useable data for this slit/order. '
                'This parameter limits the range over which the trace is fit, and may be useful if the selected slit/order '
                'would include regions without expected signal (*e.g.* bluer than the atmospheric cutoff or redder than the '
                'silicon cutoff).'
            ),
        ),
        'find_numiterfit': parset.set_parameter_definition(
            dtype=int,
            default=9,
            descr='Number of iterations to perform on the trace fitting.',
        ),
        'find_fwhm': parset.set_parameter_definition(
            dtype=[int, float],
            default=5.0,
            descr='Indicates roughly the fwhm of objects in pixels for object finding',
        ),
        'fof_link': parset.set_parameter_definition(
            dtype=[int, float],
            default=1.5,
            descr='The linking distance, in arcseconds, for the Friends of Friends algorithm to link objects across traces in Echelle spectrographs. ',
        ),
        'ech_find_max_snr': parset.set_parameter_definition(
            dtype=[int, float],
            default=1.0,
            descr=(
                'Criteria for keeping echelle objects. They must either have a maximum S/N across all the orders greater than this value '
                ' or satisfy the min_snr criteria described by the min_snr parameters. If maxnumber is set (see above) then these criteria '
                'will be applied but only the maxnumber highest (median) S/N ratio objects will be kept. '
            ),
        ),
        'ech_find_min_snr': parset.set_parameter_definition(
            dtype=[int, float],
            default=0.3,
            descr=(
                'Criteria for keeping echelle objects. They must either have a maximum S/N across all the orders greater than ech_find_max_snr,  value '
                ' or they must have S/N > ech_find_min_snr on >= ech_find_nabove_min_snr orders. If maxnumber is set (see above) then these criteria '
                'will be applied but only the maxnumber highest (median) S/N ratio objects will be kept. '
            ),
        ),
        'ech_find_nabove_min_snr': parset.set_parameter_definition(
            dtype=int,
            default=2,
            descr=(
                'Criteria for keeping echelle objects. They must either have a maximum S/N across '
                'all the orders greater than ech_find_max_snr,  value '
                ' or they must have S/N > ech_find_min_snr on >= ech_find_nabove_min_snr orders. '
                'If maxnumber is set (see above) then these criteria '
                'will be applied but only the maxnumber highest (median) S/N ratio objects will be kept.'
            ),
        ),
        'skip_second_find': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Only perform one round of object finding (mainly for quick_look)'
        ),
        'skip_final_global': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'If True, do not update initial sky to get global sky using updated noise model. This '
                'should be True for quicklook to save time. This should also be True for near-IR '
                'reductions which perform difference imaging, since there we fit sky-residuals rather '
                'than the sky itself, so there is no noise model to update. '
            ),
        ),
        'skip_skysub': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'If True, do not sky subtract when performing object finding. This should be set to '
                'True for example when running on data that is already sky-subtracted. '
                'Note that for near-IR difference imaging one still wants to remove sky-residuals via '
                'sky-subtraction, and so this is typically set to False'
            ),
        ),
        'find_negative': parset.set_parameter_definition(
            dtype=[bool],
            default=None,
            descr=(
                'Identify negative objects in object finding for spectra that are differenced. This is used to manually '
                'override the default behavior in PypeIt for object finding by setting this parameter to something other than None '
                'The default behavior is that PypeIt will search for negative object traces if background frames '
                'are present in the PypeIt file that are classified as "science" '
                '(i.e. via pypeit_setup -b, and setting bkg_id in the PypeIt file). If background frames are present '
                'that are classified as "sky", then PypeIt will NOT search for negative object traces. If one wishes '
                'to explicitly override this default behavior, set this parameter to True to find negative objects or False to ignore '
                'them.'
            ),
        ),
        'find_min_max': parset.set_parameter_definition(
            dtype=list,
            default=None,
            descr=(
                'It defines the minimum and maximum of your object in pixels in the spectral direction on the '
                'detector. It only used for object finding. This parameter is helpful if your object only '
                'has emission lines or at high redshift and the trace only shows in part of the detector.'
            ),
        ),
        'use_std_trace': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'If True, the trace of the standard star spectrum is used as a crutch for '
                'tracing the object spectra. This is useful when a direct trace is not possible '
                '(i.e., faint sources). Note that a standard star exposure must be included in your '
                'pypeit file, or the ``std_spec1d`` parameter must be set for this to work. '
                ),
        ),
        'std_spec1d': parset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=(
                'A PypeIt spec1d file of a previously reduced standard star. '
                'This can be used to trace the object spectra, but the ``use_std_trace`` '
                'parameter must be set to True. If provided, this overrides use of '
                'any standards included in your pypeit file; the standard exposures '
                'will still be reduced.'
            ),
        ),
    }

    def validate(self):
        if self.data['std_spec1d'] is not None:
            if not self.data.get('use_std_trace', True):
                raise ValueError('If you provide a standard star spectrum for tracing, you must set use_std_trace=True.')
            elif not Path(self.data['std_spec1d']).absolute().exists():
                raise ValueError(f'{self.data["std_spec1d"]} does not exist!')


class SkySubPar(parset.ParSet):
    """
    New-style parameter set for sky subtraction (replacement for SkySubPar).

    Mirrors the legacy `SkySubPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'skysub'

    parameters = {
        'bspline_spacing': parset.set_parameter_definition(
            dtype=[int, float],
            default=0.6,
            descr='Break-point spacing for the bspline sky subtraction fits.',
        ),
        'sky_sigrej': parset.set_parameter_definition(
            dtype=float,
            default=3.0,
            descr='Rejection parameter for local sky subtraction',
        ),
        'global_sky_std': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Global sky subtraction will be performed on standard stars. This should be turned '
                'off for example for near-IR reductions with narrow slits, since bright standards can '
                'fill the slit causing global sky-subtraction to fail. In these situations we go '
                'straight to local sky-subtraction since it is designed to deal with such situations'
            ),
        ),
        'no_poly': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Turn off polynomial basis (Legendre) in global sky subtraction',
        ),
        'no_local_sky': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='If True, turn off local sky model evaluation, but do fit object profile and perform optimal extraction',
        ),
        'user_regions': parset.set_parameter_definition(
            dtype=[str, list],
            default=None,
            descr=(
                'Provides a user-defined mask defining sky regions.  By '
                'default, the sky regions are identified automatically.  To '
                'specify sky regions for *all* slits, provide a comma separated '
                'list of percentages.  For example, setting user_regions = '
                ':10,35:65,80: selects the first 10%, the inner 30%, and the '
                'final 20% of *all* slits as containing sky.  Setting '
                'user_regions = user will attempt to load any SkyRegions '
                'files generated by the user via the pypeit_skysub_regions tool.'
            ),
        ),
        'mask_by_boxcar': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='In global sky evaluation, mask the sky region around the object by the boxcar radius (set in ExtractionPar).',
        ),
        'joint_fit': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Perform a simultaneous joint fit to sky regions using all available slits. '
                'Currently, this parameter is only used for IFU data reduction. Note that the '
                'current implementation does not account for variations in the instrument FWHM '
                'in different slits. This will be addressed by Issue #1660.'
            ),
        ),
        'max_mask_frac': parset.set_parameter_definition(
            dtype=float,
            default=0.80,
            descr=(
                'Maximum fraction of total pixels on a slit that can be masked by the input masks. '
                'If more than this threshold is masked the code will return zeros and throw a warning.'
            ),
        ),
        'local_maskwidth': parset.set_parameter_definition(
            dtype=float,
            default=4.0,
            descr='Initial width of the region in units of FWHM that will be used for local sky subtraction',
        ),
    }


class ExtractionPar(parset.ParSet):
    """
    New-style parameter set for extraction (replacement for ExtractionPar).

    Mirrors the legacy `ExtractionPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'extraction'

    parameters = {
        'boxcar_radius': parset.set_parameter_definition(
            dtype=[int, float],
            default=1.5,
            descr='Boxcar radius in arcseconds used for boxcar extraction',
        ),
        'skip_extraction': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Do not perform an object extraction',
        ),
        'skip_optimal': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Perform boxcar extraction only (i.e. skip Optimal and local skysub)',
        ),
        'std_prof_nsigma': parset.set_parameter_definition(
            dtype=float,
            default=30.0,
            descr='prof_nsigma parameter for Standard star extraction.  Prevents undesired rejection. '
                    'NOTE: Not consumed by the code at present.',
        ),
        'min_frac_prof': parset.set_parameter_definition(
            dtype=float,
            default=0.05,
            descr=(
                'For each spectral pixel, if the sum of the normalized object profile'
                ' across the spatial direction is less than this value,'
                ' the optimal extraction will also be masked. '
            ),
        ),
        'sn_gauss': parset.set_parameter_definition(
            dtype=[int, float],
            default=4.0,
            descr=(
                'S/N threshold for performing the more sophisticated optimal extraction which performs a '
                'b-spline fit to the object profile. For S/N < sn_gauss the code will simply optimal extract'
                'with a Gaussian with FWHM determined from the object finding.'
            ),
        ),
        'model_full_slit': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'If True local sky subtraction will be performed on the entire slit. If False, local sky subtraction will '
                'be applied to only a restricted region around each object. This should be set to True for either multislit '
                'observations using narrow slits or echelle observations with narrow slits'
            ),
        ),
        'use_2dmodel_mask': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr=(
                'Mask pixels rejected during profile fitting when extracting.'
                'Turning this off may help with bright emission lines.'
            ),
        ),
        'use_user_fwhm': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                'Boolean indicating if PypeIt should use the FWHM provided by the user '
                '(``find_fwhm`` in `FindObjPar`) for the optimal extraction. '
                'If this parameter is ``False`` (default), PypeIt estimates the FWHM for each '
                'detected object, and uses ``find_fwhm`` as initial guess.'
            ),
        ),
        'return_negative': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='If ``True`` the negative traces will be extracted and saved to disk',
        ),
    }


class Collate1DPar(parset.ParSet):
    """
    New-style parameter set for collating, coadding, and archiving 1D spectra.

    Mirrors the legacy `Collate1DPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'collate1d'
    
    valid_refframes = ['observed', 'heliocentric', 'barycentric']

    parameters = {
        'tolerance': parset.set_parameter_definition(
            dtype=[str, float, int],
            default=1.0,
            descr=(
                "The tolerance used when comparing the coordinates of objects. If two "
                "objects are within this distance from each other, they "
                "are considered the same object. If match_using is 'ra/dec' (the default) "
                "this is an angular distance. The defaults units are arcseconds but "
                "other units supported by astropy.coordinates.Angle can be used "
                "(`e.g.`, '0.003d' or '0h1m30s'). If match_using is 'pixel' this is a float."
            ),
        ),
        'dry_run': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                "If set, the script will display the matching File and Object Ids "
                "but will not flux, coadd or archive."
            ),
        ),
        'ignore_flux': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                "If set, the script will only coadd non-fluxed spectra even if flux data is present. "
                "Otherwise fluxed spectra are coadded if all spec1ds have been fluxed calibrated."
            ),
        ),
        'flux': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                "If set, the script will flux calibrate using archived sensfuncs before coadding."
            ),
        ),
        'outdir': parset.set_parameter_definition(
            dtype=str,
            default=os.getcwd(),
            descr=(
                "The path where all coadded output files and report files will be placed."
            ),
        ),
        'spec1d_outdir': parset.set_parameter_definition(
            dtype=str,
            default=None,
            descr=(
                "The path where all modified spec1d files are placed. These are only created if flux calibration or refframe correction are asked for."
            ),
        ),
        'exclude_slit_trace_bm': parset.set_parameter_definition(
            dtype=[list, str],
            default=[],
            descr=(
                "A list of slit trace bitmask bits that should be excluded."
            ),
        ),
        'exclude_serendip': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr=(
                "Whether to exclude SERENDIP objects from collating."
            ),
        ),
        'wv_rms_thresh': parset.set_parameter_definition(
            dtype=float,
            default=None,
            descr=(
                "If set, any objects with a wavelength RMS > this value are skipped, else all wavelength RMS values are accepted."
            ),
        ),
        'match_using': parset.set_parameter_definition(
            dtype=str,
            default='ra/dec',
            options=['pixel', 'ra/dec'],
            descr=(
                "Determines how 1D spectra are matched as being the same object. Must be either 'pixel' or 'ra/dec'."
            ),
        ),
        'refframe': parset.set_parameter_definition(
            dtype=str,
            default=None,
            options=valid_refframes,
            descr=(
                'Perform reference frame correction prior to coadding.  '
                f'Options are: {", ".join(valid_refframes)}'
            ),
        ),
    }


class ReducePar(parset.ParSet):
    """
    New-style parameter set for sky subtraction, object finding and extraction.

    Mirrors the legacy `ReducePar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'reduce'

    parameters = {
        'findobj': parset.set_parameter_definition(
            dtype=FindObjPar,
            default=FindObjPar(),
            descr='Parameters for the find object and tracing algorithms',
        ),
        'skysub': parset.set_parameter_definition(
            dtype=SkySubPar,
            default=SkySubPar(),
            descr='Parameters for sky subtraction algorithms',
        ),
        'extraction': parset.set_parameter_definition(
            dtype=ExtractionPar,
            default=ExtractionPar(),
            descr='Parameters for extraction algorithms',
        ),
        'slitmask': parset.set_parameter_definition(
            dtype=SlitMaskPar,
            default=SlitMaskPar(),
            descr='Parameters for slitmask',
        ),
        'cube': parset.set_parameter_definition(
            dtype=CubePar,
            default=CubePar(),
            descr='Parameters for cube generation algorithms',
        ),
        'trim_edge': parset.set_parameter_definition(
            dtype=list,
            default=[3, 3],
            descr=(
                'Trim the slit by this number of pixels left/right when performing sky subtraction'
            )
        ),
    }


class CalibrationsPar(parset.ParSet):
    """
    New-style parameter set for calibration frame groups and related settings.

    Mirrors the legacy `CalibrationsPar` in :mod:`pypeit.par.pypeitpar`.
    """

    default_key = 'calibrations'

    parameters = {
        'calib_dir': parset.set_parameter_definition(
            dtype=str,
            default='Calibrations',
            descr=(
                'The name of the directory for the processed calibration frames.  '
                'The host path for the directory is set by the redux_path (see ' 
                ':class:`~pypeit.par.pypeitpar.ReduxPar`).  Beware that success '
                'when changing the default value is not well tested!'
            ),
        ),
        'raise_chk_error': parset.set_parameter_definition(
            dtype=bool,
            default=True,
            descr='Raise an error if the calibration check fails',
        ),
        'bpm_usebias': parset.set_parameter_definition(
            dtype=bool,
            default=False,
            descr='Make a bad pixel mask from bias frames? Bias frames must be provided.',
        ),

        # TODO: Change the dict keys so that CalibrationsPar uses the FrameGroupPar.default_key?
        'biasframe': parset.set_parameter_definition(
            dtype=BiasFramePar,
            default=BiasFramePar(),
            descr='The frames and combination rules for the bias correction',
        ),
        'darkframe': parset.set_parameter_definition(
            dtype=DarkFramePar,
            default=DarkFramePar(),
            descr='The frames and combination rules for the dark-current correction',
        ),
        'scattlightframe': parset.set_parameter_definition(
            dtype=ScatteredLightFramePar,
            default=ScatteredLightFramePar(),
            descr='The frames and combination rules for the scattered light frames',
        ),
        'pixelflatframe': parset.set_parameter_definition(
            dtype=PixelFlatFramePar,
            default=PixelFlatFramePar(),
            descr='The frames and combination rules for the pixel flat',
        ),
        'illumflatframe': parset.set_parameter_definition(
            dtype=IllumFlatFramePar,
            default=IllumFlatFramePar(),
            descr='The frames and combination rules for the illumination flat',
        ),
        'lampoffflatsframe': parset.set_parameter_definition(
            dtype=LampOffFlatsFramePar,
            default=LampOffFlatsFramePar(),
            descr='The frames and combination rules for the lamp off flats',
        ),
        'slitless_pixflatframe': parset.set_parameter_definition(
            dtype=SlitlessPixFlatFramePar,
            default=SlitlessPixFlatFramePar(),
            descr='The frames and combination rules for the slitless pixel flat',
        ),
        'pinholeframe': parset.set_parameter_definition(
            dtype=PinholeFramePar,
            default=PinholeFramePar(),
            descr='The frames and combination rules for the pinholes',
        ),
        'alignframe': parset.set_parameter_definition(
            dtype=AlignFramePar,
            default=AlignFramePar(),
            descr='The frames and combination rules for the align frames',
        ),
        'arcframe': parset.set_parameter_definition(
            dtype=ArcFramePar,
            default=ArcFramePar(),
            descr='The frames and combination rules for the wavelength calibration',
        ),
        'tiltframe': parset.set_parameter_definition(
            dtype=TiltFramePar,
            default=TiltFramePar(),
            descr='The frames and combination rules for the wavelength tilts',
        ),
        'traceframe': parset.set_parameter_definition(
            dtype=TraceFramePar,
            default=TraceFramePar(),
            descr='The frames and combination rules for images used for slit tracing',
        ),
        'standardframe': parset.set_parameter_definition(
            dtype=StandardFramePar,
            default=StandardFramePar(),
            descr=(
                'The frames and combination rules for the spectrophotometric '
                'standard observations'
            ),
        ),
        'skyframe': parset.set_parameter_definition(
            dtype=SkyFramePar,
            default=SkyFramePar(),
            descr=(
                'The frames and combination rules for the sky background '
                'observations'
            ),
        ),
        'alignment': parset.set_parameter_definition(
            dtype=AlignPar,
            default=AlignPar(),
            descr='Define the procedure for the alignment of traces',
        ),
        'scattlight_pad': parset.set_parameter_definition(
            dtype=int,
            default=5,
            descr='Number of unbinned pixels to extend the slit edges by when masking the slits.',
        ),
        'flatfield': parset.set_parameter_definition(
            dtype=FlatFieldPar,
            default=FlatFieldPar(),
            descr='Parameters used to set the flat-field procedure',
        ),
        'wavelengths': parset.set_parameter_definition(
            dtype=WavelengthSolutionPar,
            default=WavelengthSolutionPar(),
            descr='Parameters used to derive the wavelength solution',
        ),
        'slitedges': parset.set_parameter_definition(
            dtype=EdgeTracePar,
            default=EdgeTracePar(),
            descr='Slit-edge tracing parameters',
        ),
        'tilts': parset.set_parameter_definition(
            dtype=WaveTiltsPar,
            default=WaveTiltsPar(),
            descr='Define how to trace the slit tilts using the trace frames',
        ),
    }


class PypeItPar(parset.ParSet):
    """
    New-style superset of parameters used by PypeIt.

    This is a single object used as a container for all the
    user-specified arguments used by PypeIt.
    """

    default_key = 'pypeit'

    parameters = {
        'rdx': parset.set_parameter_definition(
            dtype=ReduxPar,
            default=ReduxPar(),
            descr='PypeIt reduction rules.',
        ),
        'calibrations': parset.set_parameter_definition(
            dtype=CalibrationsPar,
            default=CalibrationsPar(),
            descr='Parameters for the calibration algorithms',
        ),
        'scienceframe': parset.set_parameter_definition(
            dtype=ScienceFramePar,
            default=ScienceFramePar(),
            descr='The frames and combination rules for the science observations',
        ),
        'reduce': parset.set_parameter_definition(
            dtype=ReducePar,
            default=ReducePar(),
            descr='Parameters determining sky-subtraction, object finding, and extraction',
        ),
        'flexure': parset.set_parameter_definition(
            dtype=FlexurePar,
            default=FlexurePar(),
            descr=(
                'Parameters used by the flexure-correction procedure.  Flexure '
                'corrections are not performed by default.  To turn on, either '
                "set the parameters in the 'flexure' parameter group or set 'flexure = True' "
                'in the ' "'rdx'" ' parameter group to use the default flexure-correction parameters.'
            ),
        ),
        'fluxcalib': parset.set_parameter_definition(
            dtype=FluxCalibratePar,
            default=FluxCalibratePar(),
            descr=(
                'Parameters used by the flux-calibration procedure.  Flux '
                'calibration is not performed by default.  To turn on, either '
                "set the parameters in the 'fluxcalib' parameter group or set 'fluxcalib = True' "
                'in the ' "'rdx'" ' parameter group to use the default flux-calibration parameters.'
            ),
        ),
        'coadd1d': parset.set_parameter_definition(
            dtype=Coadd1DPar,
            default=Coadd1DPar(),
            descr='Par set to control 1D coadds.  Only used in the after-burner script.',
        ),
        'coadd2d': parset.set_parameter_definition(
            dtype=Coadd2DPar,
            default=Coadd2DPar(),
            descr='Par set to control 2D coadds.  Only used in the after-burner script.',
        ),
        'sensfunc': parset.set_parameter_definition(
            dtype=SensFuncPar,
            default=SensFuncPar(),
            descr=(
                'Par set to control sensitivity function computation.  Only used in '
                'the after-burner script.'
            ),
        ),
        'telluric': parset.set_parameter_definition(
            dtype=TelluricPar,
            default=TelluricPar(),
            descr=(
                'Par set to control telluric fitting.  Only used in the '
                'pypeit_sensfunc and pypeit_telluric after-burner scripts.'
            ),
        ),
        'collate1d': parset.set_parameter_definition(
            dtype=Collate1DPar,
            default=Collate1DPar(),
            descr=(
                'Par set to control collating 1d spectra.  Only used in the '
                'after-burner script.'
            ),
        ),
    }

    @classmethod
    def from_dict(cls, cfg):
        """
        Overwrite the base class method to deal with the baseprocess
        functionality.
        """
        if 'baseprocess' not in cfg.keys():
            return super().from_dict(cfg)
        bp = cfg.pop('baseprocess')
        self = super().from_dict(cfg)
        baseproc = ProcessImagesPar.from_dict(bp)
        self.sync_processing(baseproc)
        return self

    # TODO: I'm not sure if the warning in the docstring is still valid.
    @classmethod
    def from_cfg_file(cls, cfg_file=None, merge_with=None, evaluate=True):
        """
        Construct the parameter set using a configuration file.

        Note that the following assert statement should always pass:

        .. code-block::

            default = PypeItPar()
            nofile = PypeItPar.from_cfg_file()
            assert default.data == nofile.data, 'This should always pass.'

        .. warning::

            When ``evaluate`` is true, the function runs
            :func:`~pypeit.par.util.eval_tuple` or
            :func:`~pypeit.par.util.ast_literal_eval` on all the entries in the
            `ConfigObj`_ dictionary, done using
            :func:`~pypeit.par.util.recursive_dict_evaluate`.  This has the
            potential to go haywire if the name of a parameter unintentionally
            happens to be identical to an imported or system-level function.  Of
            course, this can be useful by allowing one to define the function to
            use as a parameter, but it also means one has to be careful with the
            values that the parameters should be allowed to have.  The current
            way around this is to provide a list of strings that should be
            ignored during the evaluation, done using
            :func:`~pypeit.par.util._eval_ignore`.

        Parameters
        ----------
        cfg_file : :obj:`str`, optional
            The name of the configuration file that defines the default
            parameters.  This can be used to load a pypeit config file from a
            previous run that was constructed and output by pypeit.  This has to
            contain the full set of parameters, not just the subset you want to
            change.  For the latter, use `merge_with` to provide one or more
            config files to merge with the defaults to construct the full
            parameter set.
        merge_with : :obj:`str`, :obj:`list`, optional
            One or more config files with the modifications to either default
            parameters (`cfg_file` is None) or the parameters provided by
            `cfg_file`.  The modifications are performed in series so the list
            order of the config files is important.

        evaluate : :obj:`bool`, optional
            Evaluate the values in the config object before assigning them in
            the subsequent parameter sets.  The parameters in the config file
            are *always* read as strings, so this should almost always be true;
            however, see the warning below.

        Returns
        -------
        :class:`~pypeit.par.pypeitpar.PypeItPar`
            The instance of the parameter set.
        """
        # Get the base parameters in a ConfigObj instance
        cfg = ConfigObj(cls().to_config() if cfg_file is None else cfg_file)

        # Get the list of other configuration parameters to merge it with
        _merge_with = (
            [] if merge_with is None
            else ([merge_with] if isinstance(merge_with, str) else merge_with)
        )
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

    # TODO: This seems to be the main function that is used.  Can we consolidate
    # between `from_cfg_file` and `from_cfg_lines`?
    @classmethod
    def from_cfg_lines(cls, cfg_lines=None, merge_with=None, evaluate=True):
        """
        Construct the parameter set using the list of string lines read
        from a config file.

        Note that the following assert statement should always pass:

        .. code-block::

            default = PypeItPar()
            nofile = PypeItPar.from_cfg_lines()
            assert default.data == nofile.data, 'This should always pass.'

        .. warning::

            When ``evaluate`` is true, the function runs
            :func:`~pypeit.par.util.eval_tuple` or
            :func:`~pypeit.par.util.ast_literal_eval` on all the entries in the
            `ConfigObj`_ dictionary, done using
            :func:`~pypeit.par.util.recursive_dict_evaluate`.  This has the
            potential to go haywire if the name of a parameter unintentionally
            happens to be identical to an imported or system-level function.  Of
            course, this can be useful by allowing one to define the function to
            use as a parameter, but it also means one has to be careful with the
            values that the parameters should be allowed to have.  The current
            way around this is to provide a list of strings that should be
            ignored during the evaluation, done using
            :func:`~pypeit.par.util._eval_ignore`.

        Parameters
        ----------
        cfg_lines : :obj:`list`, optional
            A list of strings with lines read, or made to look like they are,
            from a configuration file.  This can be used to load lines from a
            previous run of pypeit that was constructed and output by pypeit.
            This has to contain the full set of parameters, not just the subset
            to change.  For the latter, leave this as the default value (None)
            and use `merge_with` to provide a set of lines to merge with the
            defaults to construct the full parameter set.
        merge_with : :obj:`tuple`, :obj:`list`, optional
            A tuple containing one more lists of strings with lines read, or
            made to look like they are, from a configuration file that should be
            merged with the lines provided by `cfg_lines`, or the default
            parameters.  The order of the lists in the tuple is important, as it
            sets the order in which the lines are merged.  Last in line has
            *highest* priority.  Or the input may be a list which will be taken
            as a single item described above.
        evaluate : :obj:`bool`, optional
            Evaluate the values in the config object before assigning them in
            the subsequent parameter sets.  The parameters in the config file
            are *always* read as strings, so this should almost always be true;
            however, see the warning below.

        Returns
        -------
        :class:`~pypeit.par.pypeitpar.PypeItPar`
            The instance of the parameter set.
        """
        # Get the base parameters in a ConfigObj instance
        cfg = ConfigObj(cls().to_config() if cfg_lines is None else cfg_lines)

        # Merge in additional parameters
        if merge_with is not None:
            # Check it is a tuple
            if isinstance(merge_with, list):
                merge_with = (merge_with,)
            if not isinstance(merge_with, tuple):
                raise PypeItError('Input merge_with must be a tuple.')
            # Proceed
            for f in merge_with:
                cfg.merge(ConfigObj(f))

        # Evaluate the strings if requested
        if evaluate:
            cfg = util.recursive_dict_evaluate(cfg)

        # Instantiate the object based on the configuration dictionary
        return cls.from_dict(cfg)
    
    def reset_all_processimages_par(self, **kwargs):
        """
        Change image processing parameter for *all* frame types.

        This function iteratively changes the value of all image processing
        parameters for all frame types in the :class:`CalibrationsPar`, as well
        as the science frames.

        Parameters
        ----------
        **kwargs:
            The list of keywords and values to change for all image processing
            parameters.

        Examples
        --------
        To turn off the slit-illumination correction for all frames:

        >>> from pypeit.spectrographs.util import load_spectrograph
        >>> spec = load_spectrograph('shane_kast_blue')
        >>> par = spec.default_pypeit_par()
        >>> par.reset_all_processimages_par(use_illumflat=False)

        """
        # Calibrations
        for _key in self['calibrations'].keys():
            if (
                isinstance(self['calibrations'][_key], parset.ParSet)
                and 'process' in self['calibrations'][_key].keys()
            ):
                for key,value in kwargs.items():
                    self['calibrations'][_key]['process'][key] = value
        # Science frame
        for _key in self.keys():
            if isinstance(self[_key], parset.ParSet) and 'process' in self[_key].keys():
                for key,value in kwargs.items():
                    self[_key]['process'][key] = value

    def sync_processing(self, proc_par):
        """
        Sync the processing of all the frame types based on the input
        :class:`~pypeit.par.pypeitpar.ProcessImagesPar` parameters.

        The parameters are merged in sequence starting from the parameter
        defaults, then including global adjustments provided by ``process``, and
        ending with the parameters that may have already been changed for each
        frame.

        This function can be used at anytime, but is most useful with the
        :class:`~pypeit.par.parset.from_dict` method where a ``baseprocess``
        group can be supplied to change the processing parameters for all frames
        away from the defaults.

        Parameters
        ----------
        proc_par : :class:`~pypeit.par.pypeitpar.ProcessImagesPar`
            Effectively a new set of default image processing parameters for all
            frames.

        Raises
        ------
        TypeError
            Raised if the provided parameter set is not an instance of
            :class:`~pypeit.par.pypeitpar.ProcessImagesPar`.
        """
        # Checks
        if not isinstance(proc_par, ProcessImagesPar):
            raise TypeError('Must provide an instance of ProcessImagesPar')

        # All the relevant ParSets are already ProcessImagesPar objects,
        # so we can work directly with the internal dictionaries.

        # Find the keys in the input that are different from the default
        default = ProcessImagesPar()
        base_diff = [key for key in proc_par.keys() if default[key] != proc_par[key]]

        # Calibration frames
        frames = [key for key in self['calibrations'].keys() if 'frame' in key]
        for f in frames:
            # Find the keys in self that are the same as the default
            frame_same = [
                key for key in proc_par.keys()
                if self['calibrations'][f]['process'].data[key] == default[key]
            ]
            to_change = list(set(base_diff) & set(frame_same))
            for key in to_change:
                self['calibrations'][f]['process'].data[key] = proc_par[key]

        # Science frames
        frame_same = [
            key for key in proc_par.keys()
            if self['scienceframe']['process'].data[key] == default[key]
        ]
        to_change = list(set(base_diff) & set(frame_same))
        for key in to_change:
            self['scienceframe']['process'].data[key] = proc_par[key]
