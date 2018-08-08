"""
Module to generate an optical model for a spectrograph.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings
import numpy

# ----------------------------------------------------------------------
# Utility functions
def conjugate_surface_transform(ray, surface_transform, forward=False):
    return numpy.matmul(surface_transform, ray) if forward else \
                numpy.matmul(surface_transform.T, ray)


def reflect(r, surface):
    """
    Reflect a (set of) ray(s) of a surface using the provided
    transformation.
    """
    _r = conjugate_surface_transform(r, surface, forward=True)
    _r[2] = -_r[2]
    return conjugate_surface_transform(r, surface)


# ----------------------------------------------------------------------
# General class for a reflection grating
class ReflectionGrating:
    """
    Doc this!
    """
    def __init__(self, ruling, tilt, roll, yaw, central_wave=None):
        self.ruling = ruling            # Grating groove spacing in lines per mm
        self.tilt = tilt                # Effective grating tilt (deg)
        self.roll = roll                # Grating roll (deg)
        self.yaw = yaw                  # Grating yaw (deg)

        self.central_wave = central_wave    # Central wavelength (angstroms)

        self.transform = self._get_grating_transform(self)

    def _get_grating_transform(self):
        """
        Taken from xidl/DEEP2/spec2d/pro/model/setup.pro ; same as in
        xidl/DEEP2/spec2d/pro/model/gsetup.pro

        Here tilt is the same as mu (!MU), roll is the same as roll3
        (!GR_YERR), and yaw is the same as o3 (!GR_ZERR).
        
        Assumes phin=0. (ie adopts the roll/yaw approach)

        MIRROR/GRATING: Better in thetax, thetay description (see above)
        for GRATING: need to add the third rotation Note the hack with
        phin + PI.  This is needed to keep the transformed x, y axes
        from flipping, wrt the original x,y.  Not a problem wrt
        reflections but if we want to work _within_ a system it is a
        problem. Cf. camera also note that this _forces_ us to work in a
        particular hemisphere, and thus we must make use of the negative
        theta as needed.
        """
        theta = -numpy.radians(self.tilt)
        rho = numpy.radians(self.roll)
        xsi = numpy.radians(self.yaw)

        cost = numpy.cos(theta)
        sint = numpy.sin(theta)
        cosr = numpy.cos(rho)
        sinr = numpy.sin(rho)
        cosx = numpy.cos(xsi)
        sinx = numpy.sin(xsi)

        return numpy.array([[ cosx*cosr,  sint*sinr+cost*sinx*cosr, -cost*sinr+sint*sinx*cosr],
                            [     -sinx,                 cost*cosx,                 sint*cosx],
                            [ cosx*sinr, -sint*cosr+cost*sinx*sinr,  cost*cosr+sint*sinx*sinr]])

    def reflect(self, r, wave=None, order=1):
        """
        Propagate an input ray for a given wavelength and order.

        wave is in angstroms
        ruling is in mm^-1

        Taken from xidl/DEEP2/spec2d/pro/model/qmodel.pro.
        """
        if wave is None and self.central_wave is None:
            raise ValueError('Must define a wavelength for the calculation.')
        if wave is None:
            warnings.warn('Using central wavelength for calculation.')
        _wave = self.central_wave if wave is None else wave
        # Transform into the grating conjugate surface
        r = conjugate_surface_transform(r, self.transform, forward=True)

        # Get the grating input angles
        alpha = -numpy.arctan(-r[1],-r[2])
        gamma = numpy.arctan(r[0],sqrt(numpy.square(r[1])+numpy.square(r[2])))

        # Use the grating equation to get the output angle
        beta = numpy.arcsin((order * 1e7*self.ruling * _wave / numpy.cos(gamma)) - numpy.sin(alpha))

        # Revert to ray vectors
        wavesign = 1-2*(_wave < 0)
        cosg = numpy.cos(gamma)
        r = numpy.array([numpy.sin(gamma),
                         numpy.sin(-beta*wavesign)*cosg,
                         numpy.cos(-beta*wavesign)*cosg ]).T

        # Return vectors transformed out of the grating conjugate surface
        return conjugate_surface_transform(r, self.transform)


# ----------------------------------------------------------------------
# Vanilla imaging spectrograph optical model
class OpticalModel:
    """
    Vanilla optical model for a imaging spectrograph.
    
    Model includes four elements:
        - Slit mask at the focal plane
        - Reflective Collimator
        - Reflective Grating
        - Refractive Camera

    Primary objective is to trace light rays from the focal plane
    position to a position in the imagine coordinate system of the
    camera.

    It is expected that each spectrograph will have its own optical
    model to perturb what is done in the vanilla model as necessary.

    .. todo:
        Provide a ParSet?

    Args:

    Attributes:

    Raises:
    """
    def __init__(self, pupil_distance, focal_r_surface, focal_r_curvature, mask_r_curvature,
                 mask_tilt_angle, mask_y_zeropoint, mask_z_zeropoint, collimator_d, collimator_r,
                 collimator_k, coll_tilt_err, coll_tilt_phi, grating, camera_tilt, camera_phi,
                 camera_focal_length, camera_distortions, imaging_rotation, optical_axis):

        # TODO: Are focal_r_surface and focal_r_curvature supposed to be
        # the same?
        self.pupil_distance = pupil_distance        # Pupil distance in mm
        self.focal_r_surface = focal_r_surface      # Radius of the image surface in mm
        self.focal_r_curvature = focal_r_curvature  # Focal-plane radius of curvature in mm

        self.mask_r_curvature = mask_r_curvature    # Mask radius of curvature in mm
        self.mask_tilt_angle = mask_tilt_angle      # Mask tilt angle in radians
        self.mask_y_zeropoint = mask_y_zeropoint    # Mask y zero point in mm
        self.mask_z_zeropoint = mask_z_zeropoint    # Mask z zero-point in mm

        self.collimator_d = collimator_d            # Collimator distance in mm
        self.collimator_r = collimator_r            # Collimator radius of curvature in mm
        self.collimator_k = collimator_k            # Collimator curvature constant
        self.coll_tilt_err = coll_tilt_err          # Collimator tilt error in radians
        self.coll_tilt_phi = coll_tilt_phi          # Collimator tilt phi in radians

        self.collimator_transform = self._collimator_transform()

        self.grating = grating                      # Grating, must be a Grating object

        self.camera_tilt = camera_tilt              # Camera angle in radians
        self.camera_phi = camera_phi                # Camera tilt phi angle in radians
        self.camera_focal_length = camera_focal_length  # Camera focal length in mm

        # Camera transmission transform is the same as a refelction
        # transform with the opposite angle: tranformation from
        # theta-x,y is better, although this is not ridiculous.
        self.camera_transform \
                = OpticalModel.get_reflection_transform(-self.camera_tilt,
                                                        self.camera_phi + numpy.pi/2)

        self.camera_distortions = camera_distortions    # Class to apply/remove camera distortions

        self.imaging_rotation = imaging_rotation    # Image coordinate system rotation in radians
        self.optical_axis = optical_axis            # Camera optical axis center (x,y) in mm

    def _collimator_transform(self):
        """
        tilt and tilt_error are in radians

        Taken from xidl/DEEP2/spec2d/pro/model/setup.pro:

        COLL ERROR (first order): In this case, we must remove the phi
        we put in, hence the more complex form
        """
        theta = 2 * self.coll_tilt_err
        phi = self.coll_tilt_phi + numpy.pi/2

        cost = numpy.cos(theta)
        sint = numpy.sin(theta)
        cosp = numpy.cos(phi)
        sinp = numpy.sin(phi)

        # TODO: Is this right?  It's what's in the IDL version...
        c = 1-(1-cost)*cosp*cosp * sinp*sinp + cost*cosp*cosp

        return numpy.array([[ 1-(1-cost)*sinp*sinp,   (1-cost)*cosp*sinp, -sinp*sint ],
                            [   (1-cost)*cosp*sinp,                    c,  cosp*sint ],
                            [            sint*sinp,           -sint*cosp,       cost ]])

    @staticmethod
    def get_reflection_transform(theta, phi):
        """
        General reflection transform.

        Taken from xidl/DEEP2/spec2d/pro/model/setup.pro:
        """
        cosp = numpy.cos(phi+numpy.pi/2)
        sinp = numpy.sin(phi+numpy.pi/2)
        cost = numpy.cos(theta)
        sint = numpy.sin(theta)
        return numpy.array([[      cosp,       sinp, numpy.zeros_like(phi)],
                            [-cost*sinp,  cost*cosp,                  sint],
                            [ sint*sinp, -sint*cosp,                  cost]])

    # TODO: Inverse operation not provided; see xidl/DEEP2/spec2d/pro/model/proj_to_mask.pro
    def project_mask_onto_plane(self, x, y, a):
        """
        Project slitmask coords (curved surface) onto a plane.
    
        Generic as long as self is properly defined
    
        Taken from xidl/DEEP2/spec2d/pro/model/mask_to_proj.pro:
    
        This is pure geometry.
    
        Args:
            x (float):
                x coordinate on mask surface in mm
            y (float):
                y coordinate on mask surface in mm
            a (float):
                position angle on curved mask surface
    
        Returns:
            float: Three floats, the x, y, and position angle of the
            coordinate in the (focal-)planar system
        """
        cosm = numpy.cos(x / self.mask_r_curvature)
        sint = numpy.sin(self.mask_tilt_angle)
        cost = numpy.cos(self.mask_tilt_angle)
        tant = sint/cost
        
        px = self.mask_r_curvature * numpy.sin(x / self.mask_r_curvature)
        py = (y - self.mask_r_curvature * tant * (1. - cosm)) * cost + self.mask_y_zeropoint
        pa = numpy.arctan((numpy.tan(numpy.radians(a)) - tant * px / self.mask_r_curvature)
                            * cost / cosm)
    
        # What follows is a small correction for the fact that the mask does
        # not lie exactly in the spherical image surface (where the
        # distortion values and gnomonic projection are calculated) and the
        # rays are arriving from the pupil image; thus, the exact locations
        # are moved _slightly_ wrt the telescope optical axis.  Note also
        # that these corrections are only calculated to first order.
        # TODO: Make this an optional correction?
        rho = numpy.sqrt(numpy.square(px) + numpy.square(py))
    
        # Spherical image and mask surface heights
        hs = self.r_surface * (1. - numpy.sqrt(1. - numpy.square(rho / self.r_surface)))
        hm = self.mask_z_zeropoint + y * sint + self.mask_r_curvature*(1. - cosm)
    
        px *= (1 - hs + hm) / self.pupil_distance
        py *= (1 - hs + hm) / self.pupil_distance
    
        return px, py, pa
    
    def telescope_plane_coo_to_unit_vector(self, x, y, a):
        """
        Convert the coordinates in the focal plane to the ray tracing unit
        vector.
    
        !! IGNORES ANGLE !!
    
        Taken from xidl/DEEP2/spec2d/pro/model/pre_grating.pro
        """
        r2 = numpy.square(x) + numpy.square(y)
        hm = self.focal_r_curvature - numpy.sqrt(numpy.square(self.focal_r_curvature)-r2)
        theta = numpy.arctan(numpy.sqrt(r2) / (self.pupil_distance - hm))
        phi = numpy.arctan(y, x)
    
        sint = numpy.sin(theta)
        return numpy.array([ numpy.cos(phi)*sint, numpy.sin(phi)*sint, numpy.cos(theta) ]).T

    def telescope_plane_coo_to_collimator_angle(self, x, y):
        """
        Convert the coordinates in the focal plane to the collimator
        angles.

        Taken from xidl/DEEP2/spec2d/pro/model/coll_angle.pro
        """

        r2 = numpy.square(x) + numpy.square(y)
        hm = self.focal_r_curvature - numpy.sqrt(numpy.square(self.focal_r_curvature)-r2)
        d2 = self.pupil_distance + self.collimator_d

        cott = (self.pupil_distance - hm) / numpy.sqrt(r2)
        k = 1. + (1. + self.collimator_k) * cott*cott
        d = d2 * (1. + self.collimator_k)

        # The following is general for conic sections.  In practice, a
        # parabola is fine.  Note the switch in sign for quadratic root
        b = numpy.sqrt(cott*cott + d2*k*(2.*self.collimator_r - d) \
                    / numpy.square(self.collimator_r - d))
        r = (self.collimator_r - d) * (b - cott) / k if self.collimator_r > d \
                    else (d - self.collimator_r) * (b + cott) / k

        # This is for a parabola:
        # r = !R_COLL * (sqrt (cott*cott + 2. * d2 / !R_COLL ) - cott)
    
        # The general conic form (important)
        r /= self.collimator_r

        return numpy.arctan(r / numpy.sqrt(1. - (1.+self.collimator_r)*r*r)), numpy.arctan(y, x)

    def mask_coo_to_grating_input_vectors(self, x, y):
        """
        Propagate rays from the mask plane to the grating.

        Taken from xidl/DEEP2/spec2d/pro/model/pre_grating.pro
        """
        if x.ndim != 1 or y.ndim != 1:
            raise ValueError('Must provide vectors.')
   
        # Project the mask coordinates onto the telescope focal plane
        px, py, pa = self.project_mask_onto_plane(x, y, numpy.zeros_like(x))

        # Project the focal-plane coordinates onto the collimator
        theta, phi = self.telescope_plane_coo_to_collimator_angle(px, py)

        # Form the reflection transform
        collimator_reflection = OpticalModel.get_reflection_transform(theta, phi)

        # Get the telescope plane unit vectors
        r = self.telescope_plane_coo_to_unit_vector(px, py)

        # Reflect off the collimator
        r = reflect(r, collimator_reflection)

        # Transform with collimator error
        return conjugate_surface_transform(r, self.collimator_transform, forward=True)

    def ics_coo_to_grating_output_vectors(self, x, y):
        """
        Revert rays from the CCD coordinates back to the grating
        output vectors.

        INPUT IS MM

        Taken from xidl/DEEP2/spec2d/pro/model/ics_post_grating.pro
        """
        # 8. Remove the mosaic transform
        cosp = numpy.cos(-self.imaging_rotation)
        sinp = numpy.sin(-self.imaging_rotation)
        _x = x - self.optical_axis[0]
        _y = y - self.optical_axis[1]
        xp =  cosp * _x + sinp * _y
        yp = -sinp * _x + cosp * _y

        # 7. Display coords are different than the x,y,z convention adopted so far:
        xp = -xp
        yp = -yp

        # 6. Convert to x,y in focal plane:
        rp = numpy.sqrt(xp*xp + yp*yp)
        theta = numpy.arctan(rp / self.camera_focal_length)
        phi = numpy.arctan(yp, xp)

        # 5. Remove the camera distortion in theta (phi unchanged)
        if self.camera_distortions is not None:
            theta = self.camera_distortions.remove_distortion(theta)

        # 4. Convert angles into unit vectors
        sint = numpy.sin(theta)
        r = numpy.array([numpy.cos(phi)*sint, numpy.sin(phi)*sint, numpy.cos(theta)]).T

        # 3. Transform through the camera system
        return conjugate_surface_transform(r, self.camera_transform)

    def grating_output_vectors_to_ics_coo(self, r):
        """
        Revert rays from the CCD coordinates back to the grating
        output vectors.

        OUTPUT IS MM

        Inverted xidl/DEEP2/spec2d/pro/model/ics_post_grating.pro
        """
        # 3. Transform through the camera system
        r = conjugate_surface_transform(r, self.camera_transform, forward=True)

        # 4. Get angles from unit vectors
        theta = numpy.arccos(r[2])
        phi = numpy.arctan(r[1]/r[0])

        # 5. Apply the camera distortion in theta (phi unchanged)
        if self.camera_distortions is not None:
            theta = self.camera_distortions.apply_distortion(theta)

        # 6. Convert angles to x,y
        tanp = numpy.tan(phi)
        x = self.camera_focal_length * numpy.tan(theta) / numpy.sqrt(1 + numpy.square(tanp))
        y = x * tanp

        # 7. Image coordinates are flipped.
        x = -x
        y = -y

        # 8. Add the mosaic transform and return the coordinates
        cosp = numpy.cos(self.imaging_rotation)
        sinp = numpy.sin(self.imaging_rotation)
        return _x*cosp - _y*sinp + self.optical_axis[0], _x*sinp + _y*cosp + self.optical_axis[1]

    def mask_to_imaging_coordinates(self, x, y, wave, order):
        """
        Convert mask coordinates in mm to detector coordinates in pixels.

        wave is in angstroms

        Taken from xidl/DEEP2/spec2d/pro/model/qmodel.pro.
        """
        # First get the grating input vectors
        r = self.mask_coo_to_grating_input_vectors(x, y)

        # Reflect the rays off the grating
        r = self.grating.reflect(r, wave, order)

        # Propagate the rays through the camera to the detector and
        # return the imaging coordinates (in mm)
        return self.grating_output_vectors_to_ics_coo(r)


# ----------------------------------------------------------------------
# General class for mapping the image plane to multiple detectors in a
# mosaic
class DetectorMap:
    """
    Base class for detector coordinate mapping.

    !! PIXEL COORDINATES ARE 1-INDEXED !!
    """
    def __init__(self):
        # Basic detector with vanilla properties

        # Number of chips
        self.nccd = 1

        # Number of pixels for each chip in each dimension
        # TODO: Currently assumes all chips are the same
        self.npix = numpy.array([2048, 2048])

        # The size of the CCD pixels in mm
        # TODO: Currently assumes all chips are the same
        self.pixel_size = 0.015

        # Nominal gap between each CCD in each dimension in mm
        # TODO: Currently assumes all chips are the same
        self.ccd_gap = numpy.array([0])

        # Width of the CCD edge in each dimension in mm
        # TODO: Currently assumes all chips are the same
        self.ccd_edge = numpy.array([1])

        # Effective size of each chip in each dimension in pixels
        # TODO: Currently assumes all chips are the same
        self.ccd_size = self.npix + (2*self.ccd_edge + self.ccd_gap)/self.pixel_size

        # Center coordinates
        origin = numpy.array([[0.,0.]])
        offset = numpy.array([[0.,0.]])
        self.ccd_center = origin * self.ccd_size[None,:] + offset
        
        # Construct the rotation matrix
        self.rotation = numpy.radians([0])
        cosa = numpy.cos(self.rotation)
        sina = numpy.sin(self.rotation)
        self.rot_matrix = numpy.array([cosa, -sina, sina, cosa]).T.reshape(self.nccd,2,2)

    def image_coordinates(self, x_pix, y_pix, detector):
        """
        Input can be an array.
        """
        if numpy.any(detector > self.nccd):
            raise ValueError('Incorrect detector number')
        
        # Reshape into vectors
        _x = numpy.atleast_1d(x_pix)
        _y = numpy.atleast_1d(y_pix)
        if _x.shape != _y.shape:
            raise ValueError('Mismatch error between x and y shape.')
        
        # Allow the detector to be coordinate specific or one value for
        # all coordinates
        _d = numpy.atleast_1d(detector-1)
        if _d.shape != _x.shape and len(_d) == 1:
            _d = numpy.repeat(_d, numpy.prod(_x.shape)).reshape(_x.shape)

        # Keep the input shape for the returned arrays
        inp_shape = None
        if len(_x.shape) > 0:
            inp_shape = _x.shape
            _d = _d.ravel()
            _x = _x.ravel()
            _y = _y.ravel()

        # Offset by the chip center
        coo = numpy.array([_x, _y]).T - self.npix[None,:]/2

        # Rotatate and offset by the CCD center
        coo = numpy.array([numpy.matmul(self.rot_matrix[d], c) for d,c in zip(_d,coo)]) \
                    + self.ccd_center[_d,:]

        # Return with the appropriate shape
        return coo[0,0] if inp_shape is None else coo[:,0].reshape(inp_shape), \
                    coo[0,1] if inp_shape is None else coo[:,1].reshape(inp_shape)

    def ccd_coordinates(self, x_img, y_img, in_mm=True):
        """
        Input can be an array.
        """
        # Reshape into vectors and convert to pixels, if necessary
        _x = numpy.atleast_1d(1e3*x_img/self.pixel_size if in_mm else x_img)
        _y = numpy.atleast_1d(1e3*y_img/self.pixel_size if in_mm else y_img)
        if _x.shape != _y.shape:
            raise ValueError('Mismatch error between x and y shape.')

        # Keep the input shape for the returned arrays
        inp_shape = None
        if len(_x.shape) > 0:
            inp_shape = _x.shape
            _x = _x.ravel()
            _y = _y.ravel()

        # Offset by the CCD center for each chip
        coo = numpy.array([_x, _y]).T[None,:,:] - self.ccd_center[:,None,:]

        # Apply the rotation matrix and offset by the chip center
        coo = numpy.array([[numpy.matmul(r.T, _c) for _c in c] \
                                for r,c in zip(self.rot_matrix, coo)]) + self.npix[None,None,:]/2

        # Determine the associated detector
        indx = numpy.all((coo > 0) & (coo <= self.npix[None,None,:]), axis=2)
        on_ndet = numpy.sum(indx, axis=0)
        if numpy.any(on_ndet == 0):
            warnings.warn('Points may not be on any detector!')
        if numpy.any(on_ndet > 1):
            warnings.warn('Points may be on more than one detector!')
        d = (numpy.arange(self.nccd)[:,None]*numpy.ones(coo.shape[1], dtype=int)[None,:])[indx]
        if d.ndim > 1:
            d = d[:,0]
        d += 1

        # Return the coordinates
        coo = coo[indx,:]
        if coo.ndim > 2:
            coo = coo[0,:,:]
        return d if inp_shape is None else d.reshape(inp_shape), \
                    coo[0,0] if inp_shape is None else coo[:,0].reshape(inp_shape), \
                    coo[0,1] if inp_shape is None else coo[:,1].reshape(inp_shape)


