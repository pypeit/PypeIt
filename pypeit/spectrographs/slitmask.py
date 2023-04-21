"""
Module to define the SlitMask class

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from collections import OrderedDict
import warnings
import numpy
from scipy import optimize, fftpack, signal
from matplotlib import pyplot

from astropy.stats import sigma_clip

from pypeit.bitmask import BitMask
from pypeit.utils import index_of_x_eq_y
from pypeit import io

from IPython import embed

class SlitMaskBitMask(BitMask):
    """
    Mask bits used for slit mask design data.
    """
    # TODO: This is overkill at the moment. Only really useful if we
    # think there may be cases where slits can have multiple purposes.
    # If they can't (and we don't think they ever will), we could
    # simplify this so that slits can only have a single type.

    # TODO: Create a script that will dynamically write the used bits
    # to a doc for readthedocs.
    def __init__(self):
        # TODO: This needs to be an OrderedDict for now to ensure that
        # the bit assigned to each key is always the same. As of python
        # 3.7, normal dict types are guaranteed to preserve insertion
        # order as part of its data model. When/if we require python
        # 3.7, we can remove this (and other) OrderedDict usage in
        # favor of just a normal dict.
        mask = OrderedDict([
                        ('ALIGN', 'Slit used for on-sky alignment'),
                      ('SCIENCE', 'Slit containts one or more targets of scientific interest.')
                           ])
        super(SlitMaskBitMask, self).__init__(list(mask.keys()), descr=list(mask.values()))


class SlitMask:
    r"""
    Generic class for a slit mask that holds the slit positions and
    IDs.

    By default no mask bits are set. Only altered if `align` or
    `science` arguments are provided.

    Args:
        corners (array-like):
            A list or numpy.ndarray with the list of coordinates. The
            object must be 2 or 3 dimensional: 1st axis is the slit,
            2nd axis are the 4 corners, 3rd axis is the x and y
            coordinate for the corner. If two dimensional, class
            assumes there is only one slit. The size of the last two
            dimensions must always be 4 and 2, respectively. The
            input order is expected to be clockwise, starting from
            the top-right corner; i.e., the order is top-right (high
            x, low y), bottom-right (low x, low y), bottom-left (low
            x, high y), top-left (high x, high y). The x coordinates
            are along the spatial direction and the y coordinates are
            long the spectral direction. The computed length
            (difference in x), width (difference in y), and position
            angle (Cartesian angle) is relative to this assumption.
        slitid (:obj:`int`, array-like, optional):
            A list of *unique* integer IDs for each slit. If None,
            just a running 0-indexed list is used. Can be a single
            integer or have shape :math:`(N_{\rm slit},)`.
        align (:obj:`bool`, `numpy.ndarray`_, optional):
            Indicates slit(s) used for on-sky alignment. Can be a
            single boolean or have shape :math:`(N_{\rm slit},)`.
        science (:obj:`bool`, `numpy.ndarray`_, optional):
            Indicates slit(s) include a target of scientific
            interest. Can be a single boolean or have shape
            :math:`(N_{\rm slit},)`.
        onsky (`numpy.ndarray`_, optional):
            1D or 2D array with on-sky metrics for each slit. The
            shape of the array must be :math:`(5,)` or :math:`(N_{\rm
            slit},5)`, and the order along rows must match slit ID
            order (see `slitid`). The five metrics per slit are (1-2)
            right ascension and declination of the slit center, (3-4)
            the slit length and width in arcseconds, and (5) the
            position angle of the slit from N through E in degrees.
        objects (`numpy.ndarray`_, optional):
            List of objects observed as a 1D or 2D array with shape
            :math:`(9,)` or :math:`(N_{\rm obj},9)`. The nine
            elements for each object is the slit id, the object ID,
            the right ascension and declination of the target, the
            object name, the object magnitude and band, and the
            object top and bottom distances from the slit edges.
            The order of the objects does not have to match that of
            the slit IDs. Also, there can be slits without objects
            and slits with multiple objects; however, objects cannot
            be provided that are not in *any* slit (i.e., the slit
            IDs in the first column of this array have to be valid).

    Attributes:
        corners (`numpy.ndarray`_):
            See above.
        id (`numpy.ndarray`_):
            See `slitid` above.
        mask (`numpy.ndarray`_):
            Mask bits selecting the type of slit.
        onsky (`numpy.ndarray`_):
            See above.
        objects (`numpy.ndarray`_):
            See above.
        slitindx (`numpy.ndarray`_):
            The index that maps from the slit data to the object
            data. For example::

                objslitdb = self.onsky[self.slitindx]

            provides the `onsky` slit data for each object with a
            shape matched to the relevant entry in
            :attr:`self.objects`.
        center (`numpy.ndarray`_):
            The geometric center of each slit.
        top (`numpy.ndarray`_):
            The top coordinate of each slit.
        bottom (`numpy.ndarray`_):
            The bottom coordinate of each slit.
        length (`numpy.ndarray`_):
            The slit length.
        width (`numpy.ndarray`_):
            The slit width.
        pa (`numpy.ndarray`_):
            The cartesian rotation angle of the slit in degrees.
        mask_radec (:obj:`tuple`, optional):
            RA, Dec (deg) of the pointing of the mask (approximate center)
        posx_pa (:obj:`float`):
            Sky PA that points to positive x (spatial) on the detector
        negx_pa (:obj:`float`):
            Sky PA that points to negative x (spatial) on the detector
        object_names (`numpy.ndarray`_):
            Object names

    Raises:
        ValueError:
            Raised if the shape of the input corners array is incorrect
            or if the number of slit IDs does not match the number of
            slits provided.
    """
    bitmask = SlitMaskBitMask()
    def __init__(self, corners, slitid=None, align=None, science=None, onsky=None, objects=None,
                 posx_pa=None, object_names=None, mask_radec=None):

        # PA
        if posx_pa is not None:
            self.posx_pa, self.negx_pa = positive_pa(posx_pa)
        else:
            self.posx_pa, self.negx_pa = None, None

        self.mask_radec = mask_radec
        self.object_names=object_names

        # TODO: Allow random input order and then fix

        # TODO: Is counter-clockwise order more natural (the order in
        # DEIMOS slitmasks is clockwise, which is why that was chosen
        # here)

        # Convert to a numpy array if it isn't already one
        _corners = numpy.asarray(corners)
        # Check the shape
        if _corners.shape[-1] != 2 or _corners.shape[-2] != 4:
            raise ValueError('Incorrect input shape.  Must provide 4 corners with x and y for '
                             'each corner.')
        # Assign corner attribute allowing for one slit on input
        # TODO: Annoyingly, numpy.atleast_3d appends a dimension at the
        # end instead of at the beginning, like what's done with
        # numpy.atleast_2d.
        self.corners = _corners.reshape(1,4,2) if _corners.ndim == 2 else _corners

        # Assign the slit IDs
        self.slitid = numpy.arange(self.corners.shape[0]) if slitid is None \
                        else numpy.atleast_1d(slitid)
        # Check the numbers match
        if self.slitid.size != self.corners.shape[0]:
            raise ValueError('Incorrect number of slit IDs provided.')
        if len(numpy.unique(self.slitid)) != len(self.slitid):
            raise ValueError('Slit IDs must be unique!')

        # Set the bitmask
        self.mask = numpy.zeros(self.nslits, dtype=self.bitmask.minimum_dtype())
        if align is not None:
            _align = numpy.atleast_1d(align)
            if _align.size != self.nslits:
                raise ValueError('Alignment flags must be provided for each slit.')
            self.mask[_align] = self.bitmask.turn_on(self.mask[_align], 'ALIGN')
        if science is not None:
            _science = numpy.atleast_1d(science)
            if _science.size != self.nslits:
                raise ValueError('Science-target flags must be provided for each slit.')
            self.mask[_science] = self.bitmask.turn_on(self.mask[_science], 'SCIENCE')

        # On-sky coordinates of the slit center
        self.onsky = None
        if onsky is not None:
            self.onsky = numpy.atleast_2d(onsky)
            if self.onsky.shape != (self.nslits,5):
                raise ValueError('Must provide sky coordinates and slit length, width, and PA '
                                 'for each slit.')

        # Expected objects in each slit
        self.objects = None
        self.slitindx = None
        if objects is not None:
            self.objects = numpy.atleast_2d(objects)
            if self.objects.shape[1] != 9:
                raise ValueError('Must provide the slit ID, object ID, sky coordinates, object name, '
                                 'object magnitude and band, top and bottom distance for each object.')
            try:
                self.slitindx = index_of_x_eq_y(self.slitid, self.objects[:,0].astype(int),
                                                strict=True)
            except:
                # Should only fault if there are slit IDs in `objects`
                # that are not in `slitid`. In that case, return a more
                # sensible error message than what index_of_x_eq_y
                # would provide.
                raise ValueError('Some slit IDs in object list not valid.')

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
                                              numpy.diff(self.corners[:,0:2,0], axis=1))).ravel()
        self.pa[self.pa < -90] += 180
        self.pa[self.pa > 90] -= 180


    def __repr__(self):
        return '<{0}: nslits={1}>'.format(self.__class__.__name__, self.nslits)

    @property
    def nslits(self):
        """The number of slits."""
        return self.slitid.size

    @property
    def alignment_slit(self):
        """Boolean array selecting the alignment slits."""
        return self.bitmask.flagged(self.mask, 'ALIGN')

    def is_alignment(self, i):
        """Check if specific slit is an alignment slit."""
        return self.bitmask.flagged(self.mask[i], 'ALIGN')

    @property
    def science_slit(self):
        """Boolean array selecting the slits with science targets."""
        return self.bitmask.flagged(self.mask, 'SCIENCE')

    def is_science(self, i):
        """Check if specific slit should have a science target."""
        return self.bitmask.flagged(self.mask[i], 'SCIENCE')


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
            provided in pixel index units (but not necessarily index
            integers).
        mask_spat (array-like):
            Slit-mask coordinates in the spatial direction, typically
            provided in mm in the focal plane.
        guess (array-like, optional):
            The initial guess for the fit parameters. Parameters are
            currently an offset (`guess[0]`) and a scale
            (`guess[1]`).
        fix (array-like, optional):
            An array of booleans indicating if the provided guess
            parameters should be fixed during the fit.
        bounds (array-like, optional):
            An array of the lower and upper bounds for the
            parameters. Must have shape :math:`(N_{\rm par},2)`, even
            if some parameters are fixed; currently :math:`N_{\rm
            par}=2`. If None, no bounds are imposed (specifically,
            bounds are set to :math:`\pm``numpy.inf`.)
        penalty (:obj:`bool`, optional):
            Include a logarithmic penalty for slits that are matched
            to multiple slits.
        maxiter (:obj:`int`, optional):
            Maximum number of fit iterations to perform. If None,
            rejection iterations are performed until no points are
            rejected. If 1, only a single fit is performed without
            any rejection; i.e., the number of rejection iterations
            is `maxiter-1`.
        maxsep (:obj:`float`, optional):
            The maximum allowed separation between the calibrated
            coordinates of the designed slit position in pixels and
            the matched trace. If None, rejection is done iteratively
            using 5-sigma clipping.
        debug (:obj:`bool`, optional):
            Show a plot of the fit residuals after each iteration.
        fit (:obj:`bool`, optional):
            Perform the fit based on the input. If False, all
            optional entries are ignored and the user has to call the
            :func:`find_best_match` function with those same entries
            to perform the fit.

    Attributes:
        trace_spat (`numpy.ndarray`_):
            Trace coordinates in the spatial direction.
        mask_spat (`numpy.ndarray`_):
            Slit-mask coordinates in the spatial direction.
        guess_par (`numpy.ndarray`_):
            The guess model parameters.
        fit_par (`numpy.ndarray`_):
            Flag that the parameters should be fit.
        bounds (tuple):
            The boundaries imposed on the fit parameters.
        par (`numpy.ndarray`_):
            The full parameter set, including any fixed parameters.
        penalty (bool):
            Flag to apply the penalty function during the fit.
        match_coo (`numpy.ndarray`_):
            The best matching coordinates of the slit mask design
            coordinates to the slit trace pixel positions.
        match_index (`numpy.ndarray`_):
            Indices in the :attr:`mask_spat` that are best matched to
            the :attr:`trace_spat` coordinates.
        match_separation (`numpy.ndarray`_):
            Signed difference between the best-matching trace and
            mask positions in trace units; negative values mean the
            best-matching mask position is larger than the trace
            position.

    Raises:
        NotImplementedError:
            Raised if the spatial positions are not 1D arrays.
    """
    def __init__(self, trace_spat, mask_spat, trace_mask=None, guess=[0.,1.], fix=[False,False],
                 bounds=None, penalty=False, maxiter=1, maxsep=None, sigma=5, debug=False,
                 fit=False):

        # Read the input coordinates and check their dimensionality
        self.trace_spat = numpy.atleast_1d(trace_spat)
        self.mask_spat = numpy.atleast_1d(mask_spat)
        if self.trace_spat.ndim > 1 or self.mask_spat.ndim > 1:
            raise NotImplementedError('SlitRegister only allows 1D arrays on input.')
        # This mask needs to be a copy so that it can be altered and
        # compared with the input
        self.trace_mask = numpy.zeros(self.trace_spat.shape, dtype=bool) if trace_mask is None \
                                else numpy.atleast_1d(trace_mask).copy()
        if self.trace_mask.shape != self.trace_spat.shape:
            raise ValueError('Trace mask must have the same shape as the spatial positions.')

        # Input and output variables instantiated during the fit
        self.guess_par = None
        self.fit_par = None
        self.bounds = None
        self.par = None
        self.penalty = None
        self.maxsep = None
        self.sigma = None
        self.match_coo = None
        self.match_index = None
        self.match_separation = None

        # Only perform the fit if requested
        if fit:
            self.find_best_match(guess=guess, fix=fix, bounds=bounds, penalty=penalty,
                                 maxiter=maxiter, maxsep=maxsep, sigma=sigma, debug=debug)

    def _setup_to_fit(self, guess, fix, bounds, penalty, maxsep, sigma):
        """Setup the necessary attributes for the fit."""
        self.guess_par = numpy.atleast_1d(guess)
        if self.guess_par.size != 2:
            raise ValueError('Must provide two guess parameters.')
        self.par = self.guess_par.copy()
        self.fit_par = numpy.invert(numpy.atleast_1d(fix))
        if self.fit_par.size != 2:
            raise ValueError('Must indicate if each of the parameters should be fixed.')
        _bounds = numpy.array([[-numpy.inf, numpy.inf], [-numpy.inf, numpy.inf]]) \
                    if bounds is None else numpy.atleast_2d(bounds)
        if _bounds.shape != (2,2):
            raise ValueError('Must provide upper and lower bounds for each parameter.')
        self.bounds = (_bounds[:,0], _bounds[:,1])
        self.penalty = penalty
        self.maxsep = maxsep
        self.sigma = sigma
                       
    def mask_to_trace_coo(self, par=None):
        """
        Compute the mask positions in the trace coordinate system.

        This is the core method that converts the mask coordinates to
        the trace coordinates; any other method that needs these data
        call this method.

        Args:
            par (array-like, optional):
                The list of (free) parameters. If any parameters are
                fixed, they should not be included, such that the
                full parameter set is::

                    self.par[self.fit_par] = par

                If None, use :attr:`par`.

        Returns:
            `numpy.ndarray`_: The expected pixel positions of the
            trace given the mask positions and the model parameters.
        """
        if par is not None:
            self.par[self.fit_par] = par
        return self.par[0] + self.mask_spat*self.par[1]
    
    def match(self, par=None, unique=False):
        """
        Match each trace to the nearest slit position based on the
        provided or internal fit parameters (using
        :func:`mask_to_trace_coo`).

        .. note::
            Even though this method returns the information
            identically meant for :attr:`match_coo`,
            :attr:`match_separation`, and :attr:`match_index`, these
            internals are *not* altered by this function.

        Args:
            par (`numpy.ndarray`_, optional):
                The parameter vector. See :func:`mask_to_trace_coo`.
            unique (:obj:`bool`, optional):
                Force the set of matched indices to be unique. This
                can lead to large separations, which can be used to
                find traced slits that are inconsistent with the mask
                design. If the function runs out of available mask
                slits to match to, it will set the remaining indices
                to -1 and the separation to 9999.

        Returns:
            Returns three `numpy.ndarray`_ objects: (1) the mask
            coordinates in the detector frame for *all* mask
            coordinates; (2) the signed minimum difference between
            each trace and any mask position; and (3) the index in
            the array of mask coordinates that provide the best match
            to each trace coordinate. The indices the latter are only
            unique if requested.
        """
        # Get the coordinates of the slit positions in the trace
        # coordinate system
        mask_pix = self.mask_to_trace_coo(par=par)
        # Calculate the separation between each trace position with every mask position
        sep = self.trace_spat[:,None] - mask_pix[None,:]
        abssep = numpy.absolute(sep)
        indx = numpy.argmin(abssep, axis=1)
        minsep = sep[numpy.arange(indx.size),indx]
        if not unique:
            # Return the minimum separation of the trace positions and
            # any slit mask position and the index of the slit mask
            # coordinate with that minimum separation
            return mask_pix, minsep, indx

        # A set of unique matches were requested
        uniq, inv, cnt = numpy.unique(indx, return_inverse=True, return_counts=True)

        while numpy.any(cnt > 1):
            # Find which coordinates have multiple matches
            multi = (cnt>1)[inv]
            multi_indx = numpy.where(multi)[0]
            # Index of the multiple match with the closest match
            _indx = numpy.argmin(numpy.atleast_2d(abssep[multi,:][:,indx[multi]]), axis=1)
            # Keep the one with the minimum match, ...
            multi[multi_indx[_indx]] = False
            # ..., find the indices of the ones that are available to
            # be matched, ...
            avail = numpy.ones(mask_pix.size, dtype=bool)
            avail[indx[numpy.invert(multi)]] = False
            # ... and find new matches for the remainder, if any
            # are available to be matched.
            if not numpy.any(avail):
                # Set dummy values for ones that could not be matched
                indx[multi] = -1
                minsep[multi] = 9999.
                break
            # Find the new, previously unmatched indices ...
            _indx = numpy.argmin(numpy.atleast_2d(abssep[multi,:][:,avail]), axis=1)
            indx[multi] = numpy.arange(mask_pix.size)[avail][_indx]
            # ... and save their separations
            minsep[multi_indx] = sep[multi_indx,indx[multi_indx]]
            # Check if the full set of indices are now unique
            uniq, inv, cnt = numpy.unique(indx, return_inverse=True, return_counts=True)

        return mask_pix, minsep, indx
    
    def minimum_separation(self, par=None):
        r"""
        Return the minimum trace and mask separation for each trace.

        This is the method that returns the residuals that are
        minimized by `scipy.optimize.least_squares`_. Those residuals
        are provided by a call to :func:`match` *without* forcing the
        matching pairs to be unique.

        The minimum separation can be penalized (if :attr:`penalty`
        is True), which multiplies each separation by :math:`2^{dN}`
        if there are non-unique matches slit-to-trace matches. Here,
        :math:`dN` is the difference between the number of traces and
        the number of uniquely matched mask positions; i.e., if two
        traces are matched to the same mask position, the separation
        is increase by a factor of 2.

        Residuals are only returned for the unmasked trace positions;
        see :attr:`trace_mask`.

        Args:
            par (`numpy.ndarray`_, optional):
                The parameter vector. See :func:`mask_to_trace_coo`.

        Returns:
            `numpy.ndarray`_: The signed difference between the trace
            coordinates and its most closely associated slit position
            (trace-mask).
        """
        # Match slits based on their separation
        mask_pix, min_sep, min_indx = self.match(par=par) #, unique=True)
        # Only return residuals for the unmasked trace positions
        gpm = numpy.invert(self.trace_mask)
        min_sep = min_sep[gpm]
        min_indx = min_indx[gpm]

        # Use the difference between the number of slits and the number
        # of unique matches to determine and apply a penalty if
        # requested
        # TODO: Multiply them all, or just the ones that are multiply
        # matched?
        if self.penalty:
            min_sep *= numpy.power(2, min_sep.size - len(numpy.unique(min_indx)))
        return min_sep

    def find_best_match(self, guess=[0.,1.], fix=[False,False], bounds=None, penalty=False,
                        maxiter=1, maxsep=None, sigma=5, debug=False):
        r"""
        Find the best match between the trace and slit-mask
        positions.

        Populates :attr:`match_coo`, :attr:`match_separation`, and
        :attr:`match_index`; the latter is also returned.

        Args:
            guess (array-like, optional):
                The initial guess for the fit parameters. Parameters
                are currently an offset (`guess[0]`) and a scale
                (`guess[1]`).
            fix (array-like, optional):
                An array of booleans indicating if the provided guess
                parameters should be fixed during the fit.
            bounds (array-like, optional):
                An array of the lower and upper bounds for the
                parameters. Must have shape :math:`(N_{\rm par},2)`,
                even if some parameters are fixed; currently
                :math:`N_{\rm par}=2`. If None, no bounds are imposed
                (specifically, bounds are set to
                :math:`\pm``numpy.inf`.)
            penalty (:obj:`bool`, optional):
                Include a logarithmic penalty for slits that are
                matched to multiple slits.
            maxiter (:obj:`int`, optional):
                Maximum number of fit iterations to perform. If None,
                rejection iterations are performed until no points
                are rejected. If 1, only a single fit is performed
                without any rejection; i.e., the number of rejection
                iterations is `maxiter-1`.
            maxsep (:obj:`float`, optional):
                The maximum allowed separation between the calibrated
                coordinates of the designed slit position in pixels
                and the matched trace. If None, rejection is done
                iteratively using sigma clipping.
            sigma (:obj:`float`, optional):
                The sigma value to use for rejection. If None, it
                will use the default set by
                `astropy.stats.sigma_clipped_stats`.
            debug (:obj:`bool`, optional):
                Show a plot of the fit residuals after each
                iteration.

        Returns:
            `numpy.ndarray`_: The index of the slit mask position
            matched to each trace position.
        """
        # Setup the parameter data for fitting and check the input
        self._setup_to_fit(guess, fix, bounds, penalty, maxsep, sigma)

        if numpy.sum(self.fit_par) == 0:
            # Nothing to fit, so just use the input parameters to
            # construct the match
            warnings.warn('No parameters to fit!')
            self.match_coo, self.match_separation, self.match_index = self.match()
            return self.match_index

        # Check the number of iterations
        if maxiter is not None and maxiter < 1:
            warnings.warn('Must perform at least one iteration; setting maxiter=1.')

        # Perform the least squares fit and save the results
        par = self.guess_par[self.fit_par]
        bounds = (self.bounds[0][self.fit_par], self.bounds[1][self.fit_par])
        result = optimize.least_squares(self.minimum_separation, par, bounds=bounds)
        if debug:
            self.show(par=result.x, maxsep=self.maxsep, sigma=self.sigma, unique=True)

        if maxiter == 1:
            # No iterations requested, so use the parameters to get the
            # matching coordinates and return the match indices
            self.match_coo, self.match_separation, self.match_index = self.match(par=result.x)
            return self.match_index

        # Rejection iterations requested
        if maxiter is None:
            # No constraint on the number of iterations
            maxiter = numpy.inf

        # Iterate the fit up to the maximum or until no points are
        # rejected
        nfit = 1
        while nfit < maxiter:
            # Find the points to reject (bad_trace)
            bad_trace = self.trace_mismatch(maxsep=self.maxsep, sigma=self.sigma)[1]
            if len(bad_trace) == 0:
                break
            # Mask the bad traces
            self.trace_mask[bad_trace] = True

            # Re-fit, starting from previous fit
            par = self.par[self.fit_par]
            result = optimize.least_squares(self.minimum_separation, par, bounds=bounds)
            if debug:
                self.show(par=result.x, maxsep=self.maxsep, sigma=self.sigma, unique=True)
            nfit += 1

        # Use the final parameter set to get the matching coordinates
        # and return the match indices
        self.match_coo, self.match_separation, self.match_index \
                = self.match(par=result.x, unique=True)
        return self.match_index

    def show(self, par=None, maxsep=None, sigma=None, unique=True, minmax=None, synced=False):
        """
        Plot the fit residuals.

        Args:
            par (`numpy.ndarray`_, optional):
                The parameter vector. See :func:`mask_to_trace_coo`.
            maxsep (:obj:`float`, optional):
                The maximum allowed separation between the calibrated
                coordinates of the designed slit position in pixels
                and the matched trace. If None, use :attr:`maxsep`;
                see :func:`find_best_match`.
            sigma (:obj:`float`, optional):
                The sigma value to use for rejection. If None, use
                :attr:`sigma`; see :func:`find_best_match`.
            unique (:obj:`bool`, optional):
                Force the set of matched indices to be unique; see
                :func:`match`.
            minmax (array-like, optional):
                A two-element array with the minimum and maximum
                coordinate value to match to the trace data; see
                :func:`trace_mismatch`.
            synced (:obj:`bool`, optional):
                The mask coordinates being matched to are synced
                left-to-right in adjacent pairs. I.e., the indices of
                left edges are all even and the indices of all right
                edges are odd.
        """
        # Set the parameters using the relevant attributes if not provided directly
        if maxsep is None:
            maxsep = self.maxsep
        if sigma is None:
            sigma = self.sigma

        # Make sure the match information is up-to-date
        self.match_coo, self.match_separation, self.match_index \
                    = self.match(par=par, unique=unique)
        # Find missing or added traces by forcing a unique match
        missing_trace, bad_trace = self.trace_mismatch(maxsep=maxsep, sigma=sigma, minmax=minmax,
                                                       synced=synced)

        # Book-keeping for good fits
        good = numpy.invert(self.trace_mask)
        good[bad_trace] = False

        # Make the plot
        pyplot.scatter(self.trace_spat[good], self.match_separation[good],
                       color='k', zorder=1, marker='.', lw=0, s=100, label='Fit')
        if numpy.any(self.trace_mask):
            pyplot.scatter(self.trace_spat[self.trace_mask],
                           self.match_separation[self.trace_mask],
                           color='0.5', zorder=1, marker='.', lw=0, s=50, label='Masked')
        if len(bad_trace) > 0:
            pyplot.scatter(self.trace_spat[bad_trace], self.match_separation[bad_trace],
                           color='C3', zorder=1, marker='x', lw=1, s=50, label='Bad Match')
        if len(missing_trace) > 0:
            pyplot.scatter(self.match_coo[missing_trace], numpy.zeros_like(missing_trace),
                           color='C1', zorder=1, marker='.', lw=0, s=100, label='Missing')
        pyplot.xlabel('Trace location (pix)')
        pyplot.ylabel('Residuals (trace-mask; pix)')
        pyplot.title('Offset = {0:.2f}; Scale = {1:.2f}; RMS = {2:.2f}'.format(
                        self.par[0], self.par[1], numpy.std(self.match_separation[good])))
        pyplot.legend()
        pyplot.show()

    def trace_mismatch(self, maxsep=None, sigma=None, minmax=None, synced=False):
        """
        Return the mismatches between the mask and trace positions.

        Based on the best-fitting (or fixed) offset and scale
        parameters, :func:`match` is executed, forcing the
        slit-mask and trace positions pairs to be uniquely matched.

        The set of slit-mask positions without a matching trace are
        identified by finding those slits in the range relevant to
        the list of trace coordinates (see `minmax`), but without a
        matching trace index.

        .. todo::
            explain synced adjustment

        The set of mask-to-trace matches are identified as "bad" if
        they meet any of the following criteria:

            - The trace has not been masked (see :attr:`trace_mask`)
            - A unique match could not be found (see :func:`match`)
            - The absolute value of the separation is larger than the
              provided `maxsep` (when `maxsep` is not None).
            - The separation is rejected by a sigma-clipping (see
              `sigma`)

        Note that there is currently no argument that disables the
        determination of bad traces. However, bad traces are simply
        returned by the method; this function changes none of the
        class attributes.

        Args:
            maxsep (:obj:`float`, optional):
                The maximum allowed separation between the calibrated
                coordinates of the designed slit position in pixels
                and the matched trace. If None, use :attr:`maxsep`;
                see :func:`find_best_match`.
            sigma (:obj:`float`, optional):
                The sigma value to use for rejection. If None, use
                :attr:`sigma`; see :func:`find_best_match`.
            minmax (array-like, optional):
                A two-element array with the minimum and maximum
                coordinate value to match to the trace data. If None,
                this is determined from :attr:`trace_spat` and the
                standard deviation of the fit residuals.
            synced (:obj:`bool`, optional):
                The mask coordinates being matched to are synced
                left-to-right in adjacent pairs. I.e., the indices of
                left edges are all even and the indices of all right
                edges are odd.
        
        Returns:
            Two `numpy.ndarray`_ objects are returned: (1) the
            indices of mask positions without a matching trace
            position and (2) the list of trace positions identified
            as "bad."
        """
        # Check parameters are available
        if self.par is None:
            raise ValueError('No parameters are available.')

        # Set the parameters using the relevant attributes if not provided directly
        if maxsep is None:
            maxsep = self.maxsep
        if sigma is None:
            sigma = self.sigma

        # Get the coordinates and the matching indices: always use the
        # existing parameters and force the matches to be unique.
        _match_coo, _match_separation, _match_index = self.match(unique=True)

        # Selection of positions included in the fit with valid match
        # indices
        gpm = numpy.invert(self.trace_mask) & (_match_index >= 0)

        # Find the overlapping region between the trace and slit-mask
        # coordinates
        if minmax is None:
            stddev = numpy.std(_match_separation[gpm])
            _minmax = numpy.array([numpy.amin(self.trace_spat) - 3*stddev,
                                  numpy.amax(self.trace_spat) + 3*stddev])
        else:
            _minmax = numpy.atleast_1d(minmax)
            if _minmax.size != 2:
                raise ValueError('`minmax` must be a two-element array.')
        overlap = (_match_coo > _minmax[0]) & (_match_coo < _minmax[1])

        # Find any large separations
        bad = self.trace_mask | (_match_index < 0)
        if maxsep is None:
            diff = numpy.ma.MaskedArray(_match_separation, mask=bad)
            kwargs = {} if sigma is None else {'sigma': sigma}
            bad = numpy.ma.getmaskarray(sigma_clip(data=diff, **kwargs))
        else:
            bad[gpm] = numpy.absolute(_match_separation[gpm]) > maxsep

        if synced:
            # Get the union of all indices with good or missing matches
            indx = numpy.array(list(set(numpy.append(numpy.where(overlap)[0],
                                                     _match_index[gpm & numpy.invert(bad)]))))
            # If the coordinate are left-right synchronized, there
            # should be an even number of indices, and the sorted
            # sequence should always have adjacent pairs the are
            # different by one
            unsynced = numpy.where(numpy.diff(indx)[::2] != 1)[0]*2
            if len(unsynced) != 0:
                offset = numpy.ones(len(unsynced), dtype=int)
                offset[indx[unsynced] % 2 == 1] = -1
                overlap[indx[unsynced]+offset] = True
            # Make sure the last index is paired
            if indx[-1] % 2 == 0:
                # Add a right
                overlap[indx[-1]+1] = True

        # Use these to construct the list of missing traces and those
        # that are unmatched to the mask
        return numpy.array(list(set(numpy.where(overlap)[0])
                             - set(_match_index[gpm & numpy.invert(bad)]))), \
                    numpy.where(bad & numpy.invert(self.trace_mask))[0]

def xc_trace(x_trace, x_design, pix_per_mm):
    """
    Use a cross-correlation to find the offset
    """
#    _x_design = -pix_per_mm*x_design[::-1]
    _x_design = pix_per_mm*x_design
    size = slit_function_length(_x_design, oversample=10)[1]
    fftsize = fftpack.next_fast_len(size)
    design_offset, design_slits_x, design_slits_y \
            = build_slit_function(_x_design, oversample=10, size=fftsize)
    trace_offset, trace_slits_x, trace_slits_y \
            = build_slit_function(x_trace, oversample=10, size=fftsize)

#    pyplot.plot(trace_slits_x, trace_slits_y)
    pyplot.plot(design_slits_x, design_slits_y)
    pyplot.scatter(_x_design, numpy.ones(_x_design.size), color='C3', marker='.', s=50, lw=0)
    pyplot.show()

    xc = signal.correlate(numpy.roll(design_slits_y, (fftsize-size)//2),
                          numpy.roll(trace_slits_y, (fftsize-size)//2), mode='full', method='fft')
    pyplot.plot(numpy.arange(xc.size), xc)
    max_xc = numpy.argmax(xc)
    pyplot.scatter(max_xc, xc[max_xc], marker='x', color='C3')
    pyplot.show()

    import pdb
    pdb.set_trace()

def slit_function_length(edges, oversample=1):
    offset = -int(numpy.floor(numpy.amin(edges)))
    return offset, (int(numpy.ceil(numpy.amax(edges)))+1+offset)*int(oversample)

def build_slit_function(edges, size=None, oversample=1, sigma=None):
    """
    Construct a unit normalized slit function
    """
    if len(edges) % 2 != 0:
        raise ValueError('Must provide synchronized set of left and right edges.')

    offset, complete_size = slit_function_length(edges, oversample=oversample)
    _size = complete_size if size is None else size

    if _size < complete_size:
        print(_size, complete_size)
        warnings.warn('Array does not cover full edge range.')

    slit_func_x = numpy.arange(_size, dtype=float)/oversample-offset
    slit_func_y = numpy.zeros(_size, dtype=float)
    for slit in edges.reshape(-1,2):
        slit_func_y[(slit_func_x > slit[0]) & (slit_func_x < slit[1])] = 1.

    return offset, slit_func_x, slit_func_y


def positive_pa(pa:float):
    """ Modify input pa to be positive (0-360)

    Args:
        pa (float): [description]

    Returns:
        [type]: [description]
    """
    # Require it be positive
    if pa < 0.:
        pa += 360.
    # Now the complement -- also require it be positive
    comp_pa = pa - 180. if pa > 180. else pa + 180.
    # Return
    return pa, comp_pa


def correct_slitpa(slitpa, maskpa):
    """ Flip 180 degree the slit PA if the value recorded
    in the slitmask design is more than +/-90 degree from the slitmask PA.

    Args:
        slitpa (:obj:`float` or `numpy.ndarray`_): position angle of the slits.
        maskpa: (:obj:`float`): position angle of the slitmask.

    Returns:
        :obj:`float` or `numpy.ndarray`_: flipped slitpa, if it is more than +/-90 from the maskpa,
        otherwise unchanged slitpa.

    """
    # Insure maskpa is positive
    maskpa, _ = positive_pa(maskpa)

    # create a new array of slitpa
    newslitpa = numpy.atleast_1d(slitpa).copy()

    for i in range(slitpa.size):
        # make slitpa positive (if it's not)
        pa, _ = positive_pa(slitpa[i])
        # flip slitpas that are >90 and <270 degrees from the maskpa
        if 90 < (pa - maskpa) < 270:
            newslitpa[i] -= 180
        # flip slitpas that are <-90 and >-270 degrees from the maskpa
        elif -90 > (pa - maskpa) > -270:
            newslitpa[i] += 180
        # check for value >= 270
        if newslitpa[i] >= 270:
            newslitpa[i] -= 360
        # check for value <= -270
        elif newslitpa[i] <= -270:
            newslitpa[i] += 360
    return newslitpa[0] if newslitpa.size == 1 else newslitpa


def load_keck_deimoslris(filename:str, instr:str):
    """ Load up the mask design info from the header
    of the file provided

    Args:
        filename (str): 
        instr (str): Name of spectrograph
            Allowed are keck_lris_xxx, keck_deimos

    Returns:
        [type]: [description]
    """
    # Open the file
    hdu = io.fits_open(filename)

    # Build the object data
    #   - Find the index of the object IDs in the slit-object
    #     mapping that match the object catalog
    mapid = hdu['SlitObjMap'].data['ObjectID']
    catid = hdu['ObjectCat'].data['ObjectID']
    indx = index_of_x_eq_y(mapid, catid)
    # .decode() assumes encoding is the default 'utf-8'
    objname = [item.strip().decode() if isinstance(item, bytes) else item.strip()
               for item in hdu['ObjectCat'].data['OBJECT']]
    # check if each objname is an ascii str (if not BinTableHDU later in the reduction will fail)
    objname = [name if name.isascii() else name.encode('ascii', 'ignore').decode() for name in objname]
    #   - Pull out the slit ID, object ID, name, object coordinates, top and bottom distance
    objects = numpy.array([hdu['SlitObjMap'].data['dSlitId'][indx].astype(int),
                        catid.astype(int),
                        hdu['ObjectCat'].data['RA_OBJ'],
                        hdu['ObjectCat'].data['DEC_OBJ'],
                        objname,
                        hdu['ObjectCat'].data['mag'],
                        hdu['ObjectCat'].data['pBand'],
                        hdu['SlitObjMap'].data['TopDist'][indx],
                        hdu['SlitObjMap'].data['BotDist'][indx]]).T
    #   - Only keep the objects that are in the slit-object mapping
    objects = objects[mapid[indx] == catid]

    # Match the slit IDs in DesiSlits to those in BluSlits
    indx = index_of_x_eq_y(hdu['DesiSlits'].data['dSlitId'], 
                            hdu['BluSlits'].data['dSlitId'],
                            strict=True)

    # PA corresponding to positive x on detector (spatial)
    posx_pa = hdu['MaskDesign'].data['PA_PNT'][-1]
    # Insure it is positive
    posx_pa, _ = positive_pa(posx_pa)

    # flip 180 degree slit PAs if are > +/-90 degrees from maskpa (posx_pa)
    slit_pas = correct_slitpa(hdu['DesiSlits'].data['slitLPA'][indx], posx_pa)

    # Instantiate the slit mask object and return it
    try:
        hdu['BluSlits'].data['slitX0']
        indices = numpy.arange(4)
    except KeyError:
        indices = numpy.arange(4)+1
    #indices = numpy.arange(4) if instr == 'keck_deimos' else numpy.arange(4)+1 
    slit_list = []
    for index in indices:
        for cdim in ['X', 'Y']:
            slit_list.append(hdu['BluSlits'].data[f'slit{cdim}{index}'])
    slitmask = SlitMask(numpy.array(slit_list).T.reshape(-1,4,2),
                                slitid=hdu['BluSlits'].data['dSlitId'],
                                align=hdu['DesiSlits'].data['slitTyp'][indx] == 'A',
                                science=hdu['DesiSlits'].data['slitTyp'][indx] == 'P',
                                onsky=numpy.array([hdu['DesiSlits'].data['slitRA'][indx],
                                                hdu['DesiSlits'].data['slitDec'][indx],
                                                hdu['DesiSlits'].data['slitLen'][indx],
                                                hdu['DesiSlits'].data['slitWid'][indx],
                                                slit_pas]).T,
                                objects=objects,
                                #object_names=hdu['ObjectCat'].data['OBJECT'],
                                mask_radec=(hdu['MaskDesign'].data['RA_PNT'][0], 
                                            hdu['MaskDesign'].data['DEC_PNT'][0]),
                                posx_pa=posx_pa)
    # Return
    return slitmask

