
**Version**: 2.0.0

=============  ================  =================  ========================================================================
Attribute      Type              Array Type         Description                                                             
=============  ================  =================  ========================================================================
``bkpt_full``  `numpy.ndarray`_  `numpy.floating`_  Full padded knot vector stored for serialization.                       
``bkpt_gpm``   `numpy.ndarray`_  `numpy.bool`_      Active-breakpoint mask; True where a breakpoint participates in the fit.
``coeff``      `numpy.ndarray`_  `numpy.floating`_  Fitted B-spline coefficients, shape (nc,).                              
``icoeff``     `numpy.ndarray`_  `numpy.floating`_  Inverse-covariance diagonal, shape (nc,).                               
``nord``       int                                  Order of the B-spline.                                                  
=============  ================  =================  ========================================================================
