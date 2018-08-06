"""
Module to define the SlitMask class
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy

class SlitMask:
    """
    Generic class for a slit mask.

    Args:
        corners (array-like):
            A list or numpy.ndarray with the list of coordinates.  The
            object must be 2 or 3 dimensional: 1st axis it the slit, 2nd
            axis are the corners, 3rd axis is the x and y coordinate for
            the corner.  If two dimensional, class assumes there is only
            one slit.  The dimensions of the last two dimensions must
            always be 4 and 2, respectively.  The input order is
            expected to be top-right, bottom-right, bottom-left,
            top-left.  The computed length, width, and position angle is
            relative to this assumption.
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
        self.length = numpy.absolute(numpy.diff(self.corners[:,0:2,0], axis=1))
        self.width = numpy.absolute(numpy.diff(self.corners[:,1:3,1], axis=1))

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


