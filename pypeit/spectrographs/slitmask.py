"""
Module to define the SlitMask class
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy
from scipy import optimize

class SlitMask:
    """
    Generic class for a slit mask that holds the slit positions and IDs.

    Args:
        corners (array-like):
            A list or numpy.ndarray with the list of coordinates.  The
            object must be 2 or 3 dimensional: 1st axis is the slit, 2nd
            axis are the 4 corners, 3rd axis is the x and y coordinate
            for the corner.  If two dimensional, class assumes there is
            only one slit.  The size of the last two dimensions must
            always be 4 and 2, respectively.  The input order is
            expected to be top-right (high x, low y), bottom-right (low
            x, low y), bottom-left (low x, high y), top-left (high x,
            high y).  The x coordinates are along the spatial direction
            and the y coordinates are long the spectral direction.  The
            computed length (difference in x), width (difference in y),
            and position angle (Cartesian angle) is relative to this
            assumption.

        slitid (array-like, optional):
            ID numbers for each slit.  If None, just a running number
            for each of the coordinates provided.

    Attributes:
        corners (numpy.ndarray):
            See above
        center (numpy.ndarray):
            The geometric center of each slit.
        top (numpy.ndarray):
            The top coordinate of each slit.
        bottom (numpy.ndarray):
            The bottom coordinate of each slit.
        length (numpy.ndarray):
            The slit length.
        width (numpy.ndarray):
            The slit width.
        pa (numpy.ndarray):
            The cartesian rotation angle of the slit in degrees.

    Raises:
        ValueError:
            Raised if the shape of the input corners array is incorrect
            or if the number of slit IDs does not match the number of
            slits provided.
    """
    def __init__(self, corners, slitid=None):

        # TODO: Allow random input order and then fix

        # Convert to a numpy array if it isn't already one
        _corners = numpy.asarray(corners)
        # Check the shape
        if _corners.shape[-1] != 2 or _corners.shape[-2] != 4:
            raise ValueError('Incorrect input shape.  Must provide 4 corners with x and y for '
                             'each corner.')
        # Assign corner attribute allowing for one slit on input
        self.corners = _corners.reshape(1,4,2) if _corners.ndim == 2 else _corners

        # Assign the slit IDs
        self.slitid = numpy.arange(self.corners.shape[0]) if slitid is None \
                        else numpy.atleast_1d(slitid)
        # Check the numbers match
        if len(self.slitid) != self.corners.shape[0]:
            raise ValueError('Incorrect number of slit IDs provided.')

        # Center coordinates
        self.center = numpy.mean(self.corners, axis=1)

        # Top and bottom (assuming the correct input order)
        self.top = numpy.mean(numpy.roll(self.corners, 1, axis=1)[:,0:2,:], axis=1)
        self.bottom = numpy.mean(self.corners[:,1:3,:], axis=1)

        # Length and width
        self.length = numpy.absolute(numpy.diff(self.corners[:,0:2,0], axis=1)).ravel()
        self.width = numpy.absolute(numpy.diff(self.corners[:,1:3,1], axis=1)).ravel()

        # Position angle
        self.pa = numpy.degrees(numpy.arctan2(numpy.diff(self.corners[:,0:2,1], axis=1),
                                numpy.diff(self.corners[:,0:2,0], axis=1)))
        self.pa[self.pa < -90] += 180
        self.pa[self.pa > 90] -= 180

    def __getitem__(self, k):
        #TODO: Add some tests to make sure this will work
        return SlitMask(self.corners[k], slitid=self.slitid[k])

    def __repr__(self):
        return '<{0}: nslits={1}>'.format(self.__class__.__name__, self.nslits)

    @property
    def nslits(self):
        return self.corners.shape[0]


class SlitRegister:
    r"""
    Match trace and slit mask positions using a linear relation.

    The class performs a least-squares fit for the following:

    .. math::
        
        X_{\rm trace} = o + s\ X_{\rm slit},

    where :math:`s` is a scale factor that nominally converts the units
    of the input slit-mask positions to pixels and :math:`o` is the
    offset in pixels.

    Assuming the trace coordinates are in pixels and the mask
    coordinates are in mm at the focal plane, the nominal scale should
    be the plate-scale of the telescope (arcsec/mm) divided by the
    plate-scale of the detector (arcsec/pixel).  The nominal offset is
    then the negative of the coordinate of the relevant detector edge
    relative to the detector center in pixels.

    Args:
        trace_spat (array-like):
            Trace coordinates in the spatial direction, typically
            provided in pixel indices.
        mask_spat (array-like):
            Slit-mask coordinates in the spatial direction, typically
            provided in mm in the focal plane.
        guess_offset (:obj:`float`, optional):
            The initial guess for the offset (:math:`o` above).  Must be
            in the same units as the trace coordinates.  If provided as
            None, the value is fixed to 0.
        guess_scale (:obj:`float`, optional):
            The initial guess for the scale (:math:`s` above).  Must
            convert the units of the mask coordiantes to the trace
            coordinates.  If provided as None, the value is fixed to 1.
        offset_limits (array-like, optional):
            An array of the lower and upper limits to allow for the
            offset during the fit.  If None, the parameter is unbounded.
        scale_limits (array-like, optional):
            An array of the lower and upper limits to allow for the
            scale during the fit.  If None, the parameter is unbounded.
        penalty (:obj:`bool`, optional):
            Include a logarithmic penalty for slits that are matched to
            multiple slits.
        fit (:obj:`bool`, optional):
            Perform the fit based on the input.  If False, all optional
            entries are ignored and the user has to call the
            :func:`find_best_match` function with those same entries to
            peform the fit.

    Attributes:
        trace_spat (numpy.ndarray):
            Trace coordinates in the spatial direction, typically
            provided in pixel indices.
        mask_spat (numpy.ndarray):
            Slit-mask coordinates in the spatial direction, typically
            provided in mm in the focal plane.
        guess_par (numpy.ndarray):
            The guess parameters, with :math:`p_0=o` and :math:`p_1=s`.
        fit_par (numpy.ndarray):
            Flag that the parameters should be fit.
        bounds (tuple):
            The boundaries imposed on the fit parameters.
        par (numpy.ndarray):
            The full parameter set, including any fixed parameters.
        penalty (bool):
            Flat to apply the penaly function during the fit.
        match_index (numpy.ndarray):
            Indices in the :attr:`mask_spat` that are best matched to
            the :attr:`trace_spat` coordiantes.
        match_separation (numpy.ndarray):
            The difference between the best-matching trace and mask
            positions in trace units.

    Raises:
        NotImplementedError:
            Raised if the spatial positions are not 1D arrays.
    """
    def __init__(self, trace_spat, mask_spat, guess_offset=None, guess_scale=None,
                 offset_limits=None, scale_limits=None, penalty=False, fit=False):

        # Read the input coordinates and check their dimensionality
        self.trace_spat = numpy.atleast_1d(trace_spat)
        self.mask_spat = numpy.atleast_1d(mask_spat)
        if self.trace_spat.ndim > 1 or self.mask_spat.ndim > 1:
            raise NotImplementedError('SlitRegister only allows 1D arrays on input.')

        # Input and output variables instantiated during the fit
        self.guess_par = None
        self.fit_par = None
        self.bounds = None
        self.par = None
        self.penalty = None
        self.match_index = None
        self.match_separation = None

        # Only perform the fit if requested
        if fit:
            self.find_best_match(guess_offset=guess_offset, guess_scale=guess_scale,
                                 offset_limits=offset_limits, scale_limits=scale_limits,
                                 penalty=penalty)
        
    def _setup_to_fit(self, guess_offset, guess_scale, offset_limits, scale_limits, penalty):
        """Setup the necessary attributes for the fit."""
        self.guess_par = numpy.array([0 if guess_offset is None else guess_offset,
                                      1 if guess_scale is None else guess_scale])
        _offset_limits = [-numpy.inf, numpy.inf] if offset_limits is None else offset_limits
        _scale_limits = [-numpy.inf, numpy.inf] if scale_limits is None else scale_limits
        self.bounds = tuple([numpy.array([o,s]) for o,s in zip(_offset_limits, _scale_limits)])
        self.penalty = penalty
        self.fit_par = numpy.array([guess_offset is not None, guess_scale is not None])
        self.par = self.guess_par.copy()
                       
    def _fill_par(self, par):
        """Replace the free parameters with the input values."""
        self.par[self.fit_par] = par
        
    def minimum_separation(self, par=None):
        r"""
        Return the minimum trace and mask separation for each trace.

        This is the function that is minimized by the optimization
        algorithm.
        
        The calculation uses the internal parameters if `par` is not
        provided.

        The minimum separation is penalized if :attr:`penalty` is True.
        The penalty multiplies each separation by :math:`2^dN` where
        :math:`dN` is the difference between the number of traces and
        the number of uniquely matched mask positions; i.e., if two
        traces are matched to the same mask position, the separation is
        increase by a factor of 2.

        Args:
            par (numpy.ndarray, optional):
                The parameter vector.

        Returns:
            numpy.ndarray: The separation in trace coordinates between
            the trace position and its most closely associated slit
            position.
        """
        # Match slits based on their separation
        min_sep, min_indx = self.match(par=par)
        # Use the difference between the number of slits and 
        # the number of unique matches to determine and apply a penalty if requested
        if self.penalty:
            min_sep *= numpy.power(2, len(self.trace_spat) - len(numpy.unique(min_indx)))
        return min_sep

    def mask_to_trace_coo(self, par=None):
        """
        Compute the mask positions in the trace coordinate system.
        """
        if par is not None:
            # Get the full parameter set including any fixed parameters
            self._fill_par(par)

        # Translate the mask coordinates to the trace coordinates
        return self.par[0] + self.mask_spat*self.par[1]
    
    def match(self, par=None):
        """
        Match each trace to the nearest slit position based on the
        provided or internal fit parameters.

        Args:
            par (numpy.ndarray, optional):
                The parameter vector.

        Returns:
            numpy.ndarray: Returns two arrays, the minimum separation
            between each trace and any slit position and the associated
            index in the slit position vector that provides that
            minimum.
        """
        mask_pix = self.mask_to_trace_coo(par=par)

        # Calculate the separation between each trace position with every mask position
        sep = numpy.absolute(self.trace_spat[:,None] - mask_pix[None,:])

        # Return the minimum separation and the match index
        return numpy.amin(sep, axis=1), numpy.argmin(sep, axis=1)
    
    def find_best_match(self, guess_offset=None, guess_scale=None, offset_limits=None,
                        scale_limits=None, penalty=False):
        r"""
        Find the best match between the trace and slit-mask positions.

        Args:
            guess_offset (:obj:`float`, optional):
                The initial guess for the offset (:math:`o` above).
                Must be in the same units as the trace coordinates.  If
                provided as None, the value is fixed to 0.
            guess_scale (:obj:`float`, optional):
                The initial guess for the scale (:math:`s` above).  Must
                convert the units of the mask coordiantes to the trace
                coordinates.  If provided as None, the value is fixed to
                1.
            offset_limits (array-like, optional):
                An array of the lower and upper limits to allow for the
                offset during the fit.  If None, the parameter is
                unbounded.
            scale_limits (array-like, optional):
                An array of the lower and upper limits to allow for the
                scale during the fit.  If None, the parameter is
                unbounded.
            penalty (:obj:`bool`, optional):
                Include a logarithmic penalty for slits that are matched
                to multiple slits.

        Returns:
            numpy.ndarray: An array of integers that match each trace to
            the provided slit positions.
        """
        # Use the input to establish how to determine the matching
        # indices
        self._setup_to_fit(guess_offset, guess_scale, offset_limits, scale_limits, penalty)

        # Check that there's something to fit
        if numpy.sum(self.fit_par) == 0:
            warnings.warn('No parameters to fit!')
            self.match_separation, self.match_index = self.match()
        else:
            # Perform the least squares fit
            par = self.guess_par[self.fit_par]
            bounds = (self.bounds[0][self.fit_par], self.bounds[1][self.fit_par])
            result = optimize.least_squares(self.minimum_separation, par, bounds=bounds)

            # Save the result
            self.match_separation, self.match_index = self.match(par=result.x)

        # Return the matching indices
        return self.match_index

    def missing_from_trace(self):
        """
        Return the indices of slits missing from the trace positions.

        Based on the best-fitting slit positions, the list of returned
        indices are those slits that should fall in the range of the
        trace coordinates provided but do not have an associated trace
        position.

        Returns:
            list: The list of missing slits.  If all slits are accounted
            for, the list will be empty.

        Raises:
            ValueError:
                Raised if the match has not been determined yet.
        """
        if self.match_index is None:
            raise ValueError('No match as been performed.  First run fit_best_match().')
        mask_pix = self.mask_to_trace_coo()
        return list(set(numpy.where((mask_pix > numpy.amin(self.trace_spat)) 
                                    & (mask_pix < numpy.amax(self.trace_spat)))[0]) 
                    - set(self.match_index))


