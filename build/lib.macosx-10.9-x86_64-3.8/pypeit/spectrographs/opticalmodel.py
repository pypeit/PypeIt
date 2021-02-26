"""
Module to generate an optical model for a spectrograph.

.. include:: ../include/links.rst

"""
import warnings
from pypeit import msgs
import numpy
import scipy
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

        self.transform = self._get_grating_transform()

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

# NOTE: Keep this around for the time-being.
    # def reflect(self, r, wave=None, order=1):
    #     """
    #     Propagate an input ray for a given wavelength and order.

    #     wave is in angstroms
    #     ruling is in mm^-1

    #     If more than one wave provided, wavelength samples are ordered
    #     along the first axis.

    #     Taken from xidl/DEEP2/spec2d/pro/model/qmodel.pro.
    #     """
    #     if wave is None and self.central_wave is None:
    #         raise ValueError('Must define a wavelength for the calculation.')
    #     if wave is None:
    #         warnings.warn('Using central wavelength for calculation.')
    #     _wave = numpy.array([self.central_wave]) if wave is None else numpy.atleast_1d(wave)
    #     if _wave.ndim > 1:
    #         raise NotImplementedError('Input wavelength must be one number or a vector.')
    #     nwave = _wave.size

    #     _r = numpy.atleast_2d(r)
    #     if _r.ndim > 2:
    #         raise NotImplementedError('Rays must be 1D for a single ray, or 2D for multiple.')

    #     # Transform into the grating conjugate surface
    #     _r = OpticalModel.conjugate_surface_transform(_r, self.transform, forward=True)

    #     # Get the grating input angles
    #     alpha = -numpy.arctan2(-_r[:,1],-_r[:,2])
    #     gamma = numpy.arctan2(_r[:,0],numpy.sqrt(numpy.square(_r[:,1])+numpy.square(_r[:,2])))

    #     # Use the grating equation to get the output angle (minus sign
    #     # in front of sin(alpha) is for reflection)
    #     beta = numpy.arcsin((order * 1e-7 * self.ruling * _wave[None,:] / numpy.cos(gamma)[:,None])
    #                         - numpy.sin(alpha)[:,None])

    #     # Revert to ray vectors
    #     wavesign = 1-2*(_wave < 0)
    #     cosg = numpy.cos(gamma)
    #     _r = numpy.array([numpy.repeat(numpy.sin(gamma), nwave).reshape(-1,nwave),
    #                       numpy.sin(-beta*wavesign[None,:])*cosg[:,None],
    #                       numpy.cos(-beta*wavesign[None,:])*cosg[:,None] ]).T

    #     if nwave == 1:
    #         # Flatten if only one wave provided
    #         _r = _r[0]

    #     # Return vectors transformed out of the grating conjugate surface
    #     return OpticalModel.conjugate_surface_transform(_r, self.transform)

    def reflect(self, r, nslits, wave=None, order=1):
        """
        Propagate an input ray for a given wavelength and order.

        wave is in angstroms
        ruling is in mm^-1

        If more than one wave provided, wavelength samples are ordered
        along the first axis.

        Taken from xidl/DEEP2/spec2d/pro/model/qmodel.pro.

        Args:
            r (numpy.ndarray):
                Rays to propagate.
            nslits (:obj:`int`):
                Number of slits
            wave (`numpy.ndarray`_):
                The wavelengths in angstroms for the propagated coordinates.
            order (:obj:`int`):
                The grating order.

        Returns:
            Rays reflected off the grating

        """
        if wave is None and self.central_wave is None:
            msgs.error('Must define a wavelength for the calculation.')
        if wave is None:
            msgs.info('Using central wavelength for calculation.')
        _wave = numpy.array([self.central_wave]) if wave is None else numpy.atleast_1d(wave)
        if _wave.ndim > 1:
            raise NotImplementedError('Input wavelength must be one number or a vector.')
        nwave = _wave.size

        _wave_arr = numpy.tile(_wave, (nslits, 1))

        _r = numpy.atleast_2d(r)
        if _r.ndim > 2:
            raise NotImplementedError('Rays must be 1D for a single ray, or 2D for multiple.')

        # Transform into the grating conjugate surface
        _r = OpticalModel.conjugate_surface_transform(_r, self.transform, forward=True)

        # Get the grating input angles
        alpha = -numpy.arctan2(-_r[:, 1], -_r[:, 2])
        gamma = numpy.arctan2(_r[:, 0], numpy.sqrt(numpy.square(_r[:, 1]) + numpy.square(_r[:, 2])))

        # Use the grating equation to get the output angle (minus sign
        # in front of sin(alpha) is for reflection)
        beta = (numpy.arcsin((order * 1e-7 * self.ruling * _wave_arr.ravel() / numpy.cos(gamma))
                             - numpy.sin(alpha))).reshape(_wave_arr.shape).T

        # Revert to ray vectors
        wavesign = 1 - 2 * (_wave_arr.T < 0)
        _r = numpy.array([numpy.sin(gamma),
                          numpy.sin(-beta * wavesign).T.flatten() * numpy.cos(gamma),
                          numpy.cos(-beta * wavesign).T.flatten() * numpy.cos(gamma)]).T

        #if nwave == 1:
        #    # Flatten if only one wave provided
        #    _r = _r[0]

        # Return vectors transformed out of the grating conjugate surface
        return OpticalModel.conjugate_surface_transform(_r, self.transform)


# ----------------------------------------------------------------------
# Vanilla imaging spectrograph optical model
class OpticalModel:
    """
    Vanilla optical model for an imaging spectrograph.
    
    Model includes four elements:
        - Slit mask at the focal plane
        - Reflective Collimator
        - Reflective Grating
        - Refractive Camera

    Primary objective is to trace light rays from the focal plane
    position to a position in the imagine coordinate system of the
    camera.  See :func:`mask_to_imaging_coordinates`.

    It is expected that each spectrograph will have its own optical
    model to perturb what is done in the vanilla model as necessary.

    .. todo:
        - Provide a ParSet for the arguments of the function?
        - I say this is vanilla, but I'm not sure how much of this is
          still DEIMOS specific.

    Args:
        pupil_distance (float):
            Pupil distance in mm
        focal_r_surface (float):
            Radius of the image surface in mm
        focal_r_curvature (float):
            Focal-plane radius of curvature in mm
        mask_r_curvature (float):
            Mask radius of curvature in mm
        mask_tilt_angle (float):
            Mask tilt angle in radians
        mask_y_zeropoint (float):
            Mask y zero point in mm
        mask_z_zeropoint (float):
            Mask z zero point in mm
        collimator_d (float):
            Collimator distance in mm
        collimator_r (float):
            Collimator radius of curvature in mm
        collimator_k (float):
            Collimator curvature constant
        coll_tilt_err (float):
            Collimator tilt error in radians
        coll_tilt_phi (float):
            Collimator tilt phi in radians
        grating (:class:`ReflectionGrating`):
            Grating object that evaluates the grating equation and
            performs the grating reflection.
        camera_tilt (float):
            Camera angle in radians
        camera_phi (float):
            Camera tilt phi angle in radians
        camera_focal_length (float):
            Camera focal length in mm
        camera_distortions (object):
            Class to apply/remove camera distortions.  Can be None.  If
            provided, must have `remove_distortion` and
            `apply_distortion` methods.
        imaging_rotation (float):
            Image coordinate system rotation in radians
        optical_axis (numpy.ndarray):
            Camera optical axis center (x,y) in mm

    Attributes:
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
                                                        self.camera_phi + numpy.pi)

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
    def conjugate_surface_transform(ray, surface_transform, forward=False):
        """
        Transform a ray by a surface.

        Taken from xidl/DEEP2/spec2d/pro/model/gen_xfm.pro

        .. todo::
            I'm trying to mimic what was done in the IDL code,
            meaning the matrix ordering may not make sense...

        Args:
            ray (numpy.ndarray):
                The rays to tranform.  If more than one ray is provided,
                they must be organized along the last axis.  I.e., for
                two rays, the shape of ray should be (2,3).
            surface_transform (numpy.ndarray):
                The transform for the surface.  For example, see
                :func:`OpticalModel.get_reflection_transform`.  For more
                than one surface, the surface matrix must be organized
                along the last two axes.  I.e., for two surfaces, the
                shape of surface should be (2,3,3).
        Returns:
            numpy.ndarray: The array with the reflected arrays.
        """
        return numpy.einsum('...ji,...i', surface_transform, ray) if forward else \
                    numpy.einsum('...ij,...i', surface_transform, ray)

    @staticmethod
    def reflect(r, surface):
        """
        Reflect a (set of) ray(s) of a surface using the provided
        transformation.

        .. todo::
            (KBW) I'm trying to mimic what was done in the IDL code,
            meaning the matrix ordering may not make sense...

        Args:
            r (numpy.ndarray):
                The rays to reflect of the defined surface.  If more
                than one ray is provided, they must be organized along
                the last axis.  I.e., for two rays, the shape of r
                should be (2,3).
            surface (numpy.ndarray):
                The transform for the surface used in the reflection.
                For example, see
                :func:`OpticalModel.get_reflection_transform`.  For more
                than one surface, the surface matrix must be organized
                along the last two axes.  I.e., for two surfaces, the
                shape of surface should be (2,3,3).

        Returns:
            numpy.ndarray: The array with the reflected arrays.
        """
        _r = OpticalModel.conjugate_surface_transform(r, surface, forward=True)
        _r[...,2] *= -1
        return OpticalModel.conjugate_surface_transform(_r, surface)

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
        transform = numpy.array([[      cosp,       sinp, numpy.zeros_like(phi)],
                                 [-cost*sinp,  cost*cosp,                  sint],
                                 [ sint*sinp, -sint*cosp,                  cost]])

        # If theta and phi are vectors move the axes so that, e.g.,
        # transform[0] is the transform for (theta[0],phi[0])
        return numpy.moveaxis(transform, -1, 0) if transform.ndim > 2 else transform

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
        hs = self.focal_r_surface * (1. - numpy.sqrt(1. - numpy.square(rho/self.focal_r_surface)))
        hm = self.mask_z_zeropoint + y * sint + self.mask_r_curvature*(1. - cosm)
    
        px *= (1 - (hs - hm) / self.pupil_distance)
        py *= (1 - (hs - hm) / self.pupil_distance)

        return px, py, pa
    
    def telescope_plane_coo_to_unit_vector(self, x, y): #, a):
        """
        Convert the coordinates in the focal plane to the ray tracing unit
        vector.
    
        !! IGNORES ANGLE !!
    
        Taken from xidl/DEEP2/spec2d/pro/model/pre_grating.pro
        """
        r2 = numpy.square(x) + numpy.square(y)
        hm = self.focal_r_curvature - numpy.sqrt(numpy.square(self.focal_r_curvature)-r2)
        theta = numpy.arctan(numpy.sqrt(r2) / (self.pupil_distance - hm))
        phi = numpy.arctan2(y, x)
    
        sint = numpy.sin(theta)
        return numpy.array([ numpy.cos(phi)*sint, numpy.sin(phi)*sint, numpy.cos(theta) ]).T

    def telescope_plane_coo_to_collimator_angle(self, x, y):
        """
        Convert the coordinates in the focal plane to the collimator
        angles.

        Taken from xidl/DEEP2/spec2d/pro/model/coll_angle.pro
        """
        r2 = numpy.square(x) + numpy.square(y)
        fr2 = numpy.square(self.focal_r_curvature)
        hm = self.focal_r_curvature - numpy.sqrt(fr2-r2)
        d2 = self.pupil_distance + self.collimator_d

        cott = (self.pupil_distance - hm) / numpy.sqrt(r2)
        cott2 = numpy.square(cott)
        k = 1. + (1. + self.collimator_k) * cott2
        d = d2 * (1. + self.collimator_k)

        # The following is general for conic sections.  In practice, a
        # parabola is fine.  Note the switch in sign for quadratic root
        b = numpy.sqrt(cott2 + d2*k*(2.*self.collimator_r - d) \
                    / numpy.square(self.collimator_r - d))
        r = (self.collimator_r - d) * (b - cott) / k if self.collimator_r > d \
                    else (d - self.collimator_r) * (b + cott) / k

        # This is for a parabola:
        # r = !R_COLL * (sqrt (cott*cott + 2. * d2 / !R_COLL ) - cott)
    
        # The general conic form (important)
        r /= self.collimator_r

        return numpy.arctan(r / numpy.sqrt(1. - (1.+self.collimator_k)*r*r)), numpy.arctan2(y, x)

    def mask_coo_to_grating_input_vectors(self, x, y):
        """
        Propagate rays from the mask plane to the grating.

        Taken from xidl/DEEP2/spec2d/pro/model/pre_grating.pro
        """
        _x = numpy.atleast_1d(x)
        _y = numpy.atleast_1d(y)

        # Project the mask coordinates onto the telescope focal plane
        px, py, pa = self.project_mask_onto_plane(_x, _y, numpy.zeros_like(_x))

        # Project the focal-plane coordinates onto the collimator
        theta, phi = self.telescope_plane_coo_to_collimator_angle(px, py)

        # Form the reflection transform (a1)
        collimator_reflection = OpticalModel.get_reflection_transform(theta, phi)

        # Get the telescope plane unit vectors
        r = self.telescope_plane_coo_to_unit_vector(px, py)

        # Reflect off the collimator
        r = OpticalModel.reflect(r, collimator_reflection)

        # Transform with collimator error (e1)
        return OpticalModel.conjugate_surface_transform(r, self.collimator_transform,
                                                        forward=True)

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
        xp *= -1
        yp *= -1

        # 6. Convert to x,y in focal plane:
        rp = numpy.sqrt(xp*xp + yp*yp)
        theta = numpy.arctan(rp / self.camera_focal_length)
        phi = numpy.arctan2(yp, xp)

        # 5. Remove the camera distortion in theta (phi unchanged)
        if self.camera_distortions is not None:
            theta = self.camera_distortions.remove_distortion(theta)

        # 4. Convert angles into unit vectors
        sint = numpy.sin(theta)
        r = numpy.array([numpy.cos(phi)*sint, numpy.sin(phi)*sint, numpy.cos(theta)]).T

        # 3. Transform through the camera system
        return OpticalModel.conjugate_surface_transform(r, self.camera_transform)

    def grating_output_vectors_to_ics_coo(self, r, sign=1):
        """
        Revert rays from the CCD coordinates back to the grating
        output vectors.

        There's a sign degeneracy going this way, so it must be defined.

        OUTPUT IS MM

        Inverted xidl/DEEP2/spec2d/pro/model/ics_post_grating.pro
        """
        # 3. Transform through the camera system.
        _r = OpticalModel.conjugate_surface_transform(r, self.camera_transform, forward=True)

        # 4. Get angles from unit vectors
        theta = numpy.arccos(_r[...,2])
        phi = numpy.arctan2(_r[...,1], _r[...,0])

        # 5. Apply the camera distortion in theta (phi unchanged)
        if self.camera_distortions is not None:
            theta = self.camera_distortions.apply_distortion(theta)

        # 6. Convert angles to x,y
        # WARNING: There's a +/- 90 deg limit on imaging_rotation when
        # this will work...
        tanp = numpy.tan(phi)
        xp = self.camera_focal_length * numpy.tan(theta) / numpy.sqrt(1 + numpy.square(tanp))
        yp = xp * tanp

        # 7. Image coordinates are flipped.
        xp *= sign
        yp *= sign

        # 8. Add the mosaic transform and return the coordinates
        cosp = numpy.cos(-self.imaging_rotation)
        sinp = numpy.sin(-self.imaging_rotation)
        return xp*cosp - yp*sinp + self.optical_axis[0], xp*sinp + yp*cosp + self.optical_axis[1]

    def pre_grating_vectors(self, x, y, amap, npoints=1):
        """
        Propagate rays from the mask plane to the grating, by interpolating a pre-grating
        map (amap).

        This should replace :attr:`mask_coo_to_grating_input_vectors`

        Taken from DEEP2/spec2d/pro/model/qmodel.pro

        Args:
            x (`numpy.ndarray`_):
                The x coordinates in the slit mask in mm.
            y (`numpy.ndarray`_):
                The y coordinates in the slit mask in mm.
            amap (`FITS_rec`):
                pre-grating map
            npoints (:obj:`int`):
                Size of the spectral direction

        Returns:
            `numpy.ndarray`_: Rays propagated from mask plane to grating.

        """

        xmm = numpy.tile(x, (npoints, 1))
        ymm = numpy.tile(y, (npoints, 1))

        sx = amap['tanx'].squeeze().T.shape[0]
        sy = amap['tanx'].squeeze().T.shape[1]

        # create a set of indices to interpolate amap into
        xindx = (xmm - amap['xmin']) / amap['xstep']
        yindx = (ymm - amap['ymin']) / amap['ystep']

        # preparing input and output coordinates for interpolation
        _x, _y = numpy.meshgrid(amap['xarr'].squeeze(), amap['yarr'].squeeze())
        out_coo = numpy.column_stack((xmm.ravel(), ymm.ravel()))
        in_coo = numpy.column_stack((_x.ravel(), _y.ravel()))

        interp = scipy.interpolate.CloughTocher2DInterpolator(in_coo, amap['tanx'].ravel(),
                                                              fill_value=-1e10)
        tanxx = interp(out_coo).reshape(xindx.shape)

        interp = scipy.interpolate.CloughTocher2DInterpolator(in_coo, amap['tany'].ravel(),
                                                              fill_value=-1e10)
        tanyy = interp(out_coo).reshape(xindx.shape)

        whbad = (xindx < 4) | (xindx > (sx - 4)) | (yindx < 4) | (yindx > (sy - 4))
        tanxx[whbad] = -1e10
        tanyy[whbad] = -1e10

        rr_2 = (-1. / numpy.sqrt(1. + numpy.square(tanxx).T + numpy.square(tanyy).T)).ravel()
        rr = numpy.array([rr_2 * tanxx.T.ravel(), rr_2 * tanyy.T.ravel(), rr_2]).T

        return rr

    def post_grating_vectors_to_ics_coo(self, r, bmap, nslits, npoints):
        """
        Revert rays from post-grating output vectors to CCD coordinates, by interpolating a post-grating
        map (bmap).

        This should replace :attr:`grating_output_vectors_to_ics_coo`

        Taken from DEEP2/spec2d/pro/model/qmodel.pro

        Args:
            r (`numpy.ndarray`_):
                Rays to be transformed
            bmap (`FITS_rec`):
                post-grating map
            nslits (:obj:`int`):
                Number of slits
            npoints (:obj:`int`):
                Size of the spectral direction

        Returns:
            Two `numpy.ndarray`_:  Detector image plane coordinates in pixels

        """

        tanxi = r[:, 0] / r[:, 2]
        tanyi = r[:, 1] / r[:, 2]

        # use grids to determine values of xics, yics from tanxi, tanyi
        xindx = (tanxi - bmap['txmin']) / bmap['txstep']
        yindx = (tanyi - bmap['tymin']) / bmap['tystep']

        sx = bmap['gridx'].squeeze().T.shape[0]
        sy = bmap['gridx'].squeeze().T.shape[1]

        # preparing input and output coordinates for interpolation
        indx = numpy.logical_not(bmap['gridx'].squeeze() == -1e10)
        _x, _y = numpy.meshgrid(numpy.arange(sx), numpy.arange(sy))
        out_coo = numpy.column_stack((xindx.ravel(), yindx.ravel()))
        in_coo = numpy.column_stack((_x.ravel()[indx.ravel()], _y.ravel()[indx.ravel()]))

        interp = scipy.interpolate.CloughTocher2DInterpolator(in_coo, bmap['gridx'].ravel()[indx.ravel()],
                                                              fill_value=-1e10)
        xics = interp(out_coo)

        interp = scipy.interpolate.CloughTocher2DInterpolator(in_coo, bmap['gridy'].ravel()[indx.ravel()],
                                                              fill_value=-1e10)
        yics = interp(out_coo)

        whbad = (xindx < 4) | (xindx > (sx - 4)) | (yindx < 4) | (yindx > (sy - 4))
        xics[whbad] = -1e10
        yics[whbad] = -1e10
        # this condition may not be necessary. We should not have values > 1e9
        wh = (numpy.abs(xics) > 1e9) | (numpy.abs(yics) > 1e9)
        xics[wh] = -1e4
        yics[wh] = -1e4

        if npoints > 1:
            xics = xics.reshape(nslits, npoints)
            yics = yics.reshape(nslits, npoints)

        return xics, yics

    def mask_to_imaging_coordinates(self, x, y, amap, bmap, nslits, wave, order):
        """
        Convert mask coordinates in mm to detector coordinates in pixels.

        wave is in angstroms

        If more than one wavelength is provided, wavelength samples are
        ordered along the first axis.

        Taken from xidl/DEEP2/spec2d/pro/model/qmodel.pro.

        Args:
            x (`numpy.ndarray`_):
                The x coordinates in the slit mask in mm.
            y (`numpy.ndarray`_):
                The y coordinates in the slit mask in mm.
            amap (`FITS_rec`):
                pre-grating map
            bmap (`FITS_rec`):
                post-grating map
            nslits (:obj:`int`):
                Number of slits
            wave (`numpy.ndarray`_):
                The wavelengths in angstroms for the propagated coordinates.
            order (:obj:`int`):
                The grating order.

        Returns:
            Two `numpy.ndarray`_: Detector image plane coordinates in pixels
        """

        npoints = 1 if wave is None else numpy.atleast_1d(wave).shape[0]

        # First get the grating input vectors
        # r = self.mask_coo_to_grating_input_vectors(x, y)
        r = self.pre_grating_vectors(x, y, amap, npoints=npoints)

        # Reflect the rays off the grating
        r = self.grating.reflect(r, nslits, wave, order)

        # Propagate the rays through the camera to the detector and
        # return the imaging coordinates (in mm)
        # return self.grating_output_vectors_to_ics_coo(r, sign=1 - 2 * (x < 0))
        return self.post_grating_vectors_to_ics_coo(r, bmap, nslits, npoints)

# NOTE: Keep this around for the time-being.
    # def mask_to_imaging_coordinates(self, x, y, wave, order):
    #     """
    #     Convert mask coordinates in mm to detector coordinates in pixels.

    #     wave is in angstroms

    #     If more than one wavelength is provided, wavelength samples are
    #     ordered along the first axis.

    #     Taken from xidl/DEEP2/spec2d/pro/model/qmodel.pro.
    #     """
    #     # First get the grating input vectors
    #     r = self.mask_coo_to_grating_input_vectors(x, y)

    #     # Reflect the rays off the grating
    #     r = self.grating.reflect(r, wave, order)

    #     # Propagate the rays through the camera to the detector and
    #     # return the imaging coordinates (in mm)
    #     return self.grating_output_vectors_to_ics_coo(r, sign=1-2*(x<0))


# ----------------------------------------------------------------------
# Detector coordinates
class DetectorMap:
    """
    General class for mapping the image plane to the pixel coordinates
    for multiple detectors in a mosaic.

    !! PIXEL COORDINATES ARE 1-INDEXED !!

    All CCDs in the detector are assumed to have the same size.

    .. todo:
        - Allow for instantiation arguments.
        - Remove ccd_gap and ccd_edge from attributes?

    Attributes:
        nccd (int):
            Number of CCD chips in the detector.
        npix (numpy.ndarray):
            Number of pixels in each dimension (x,y) of all CCDs.
        pixel_size (float):
            Pixel size in mm.
        ccd_gap (numpy.ndarray):
            The nominal gap between each chip in (x,y) in mm.
        ccd_edge (numpy.ndarray):
            The width of the CCD edge in (x,y) in mm.
        ccd_size (numpy.ndarray):
            The nominal size in number of pixels of each CCD accounting
            for gap, edge width, and number of pixels.
        ccd_center (numpy.ndarray):
            The (x,y) center of each CCD in number of pixels, accounting
            for the per CCD offset.
        rotation (numpy.ndarray):
            The rotation of each CCD.
        rot_matrix (numpy.ndarray):
            The rotation matrix used for each CCD.
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
        self.ccd_gap = numpy.array([0, 0])

        # Width of the CCD edge in each dimension in mm
        # TODO: Currently assumes all chips are the same
        self.ccd_edge = numpy.array([1, 1])

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

    def image_coordinates(self, x_pix, y_pix, detector=1, in_mm=True):
        """
        Convert the provided (1-indexed) pixel coordinates into
        coordinates in the image plane.

        Args:
            x_pix (:obj:`float` or array-like):
                Pixel coordinate in x (1-indexed) on the detector
            y_pix (:obj:`float` or array-like):
                Pixel coordinate in y (1-indexed) on the detector
            detector (:obj:`int` or array-like, optional):
                Relevant detector for the pixel coordinates.  Default is
                1.  Can be a single detector used for all coordinates or
                an array that provides the detector for each (x,y)
                coordinate.
            in_mm (:obj:`bool`, optional):
                Return the coordinates in mm.

        Returns:
            float, numpy.ndarray: Returns two objects with the x and y
            coordinates.  Return object type is based on the input.

        Raises:
            ValueError:
                Raised if the detector number is not valid or if the
                input x and y arrays do not have the same shape.
        """
        # Reshape into vectors
        _x = numpy.atleast_1d(x_pix)
        _y = numpy.atleast_1d(y_pix)
        if _x.shape != _y.shape:
            raise ValueError('Mismatch error between x and y shape.')
        
        # Allow the detector to be coordinate specific or one value for
        # all coordinates
        _d = numpy.atleast_1d(detector)-1
        if _d.shape != _x.shape and len(_d) == 1:
            _d = numpy.full(_x.shape, detector-1, dtype=int)
        if numpy.any((_d >= self.nccd) | (_d < 0)):
            raise ValueError('Incorrect detector number')

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

        x_img = coo[0,0] if inp_shape is None else coo[:,0].reshape(inp_shape)
        y_img = coo[0,1] if inp_shape is None else coo[:,1].reshape(inp_shape)

        # Return with the appropriate shape
        return x_img * self.pixel_size if in_mm else x_img, \
                    y_img * self.pixel_size if in_mm else y_img

    def ccd_coordinates(self, x_img, y_img, in_mm=True):
        """
        Convert the provided coordinates in the image plane to
        (1-indexed) pixel coordinates on the relevant detector.

        Args:
            x_img (:obj:`float` or array-like):
                Image coordinates in x
            y_img (:obj:`float` or array-like):
                Image coordinates in y
            in_mm (:obj:`bool`, optional):
                The input coordinates are provided in mm, not pixels.

        Returns:
            float, numpy.ndarray: Returns three objects, the detector
            associated with each coordinates and the x and y pixel
            coordinates on that detector.  Return object type (float or
            array) is based on the input.

        Raises:
            ValueError:
                Raised if the input x and y arrays do not have the same
                shape.
        """

        # Reshape into vectors and convert to pixels, if necessary
        _x = numpy.atleast_1d(x_img/self.pixel_size if in_mm else x_img)
        _y = numpy.atleast_1d(y_img/self.pixel_size if in_mm else y_img)
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

        # Determine the associated detector (1-indexed)
        indx = numpy.all((coo > 0) & (coo <= self.npix[None,None,:]), axis=2)
        on_ndet = numpy.sum(indx, axis=0)
        if numpy.any(on_ndet == 0):
            warnings.warn('Points may not be on any detector!')
        if numpy.any(on_ndet > 1):
            warnings.warn('Points may be on more than one detector!')
        d = numpy.amax(numpy.arange(self.nccd)[:,None]*indx, axis=0) + 1
        d[numpy.sum(indx, axis=0) == 0] = -1

        # Pull out the coordinates for the correct detector
        coo = numpy.array([coo[_d-1,i,:] if _d > 0 else numpy.full(coo.shape[2:],-1.) 
                                for i,_d in enumerate(d)])

        # Return the coordinates
        return d if inp_shape is None else d.reshape(inp_shape), \
                    coo[0,0] if inp_shape is None else coo[:,0].reshape(inp_shape), \
                    coo[0,1] if inp_shape is None else coo[:,1].reshape(inp_shape)


