
**Version**: 2.0.0

=============  ================  =================  =======================================================================================================================
Attribute      Type              Array Type         Description                                                                                                            
=============  ================  =================  =======================================================================================================================
``basis``      `numpy.ndarray`_  `numpy.floating`_  Polynomial basis matrix, shape (N, npoly); stored only when fit was called with a numpy array basis (funcname is None).
``bkpt_full``  `numpy.ndarray`_  `numpy.floating`_  Full padded knot vector stored for serialization.                                                                      
``bkpt_gpm``   `numpy.ndarray`_  `numpy.bool`_      Active-breakpoint mask; True where a breakpoint participates in the fit.                                               
``coeff``      `numpy.ndarray`_  `numpy.floating`_  Fitted 2D B-spline coefficients, shape (nc, npoly).                                                                    
``funcname``   str                                  Polynomial family name (e.g. ``legendre``); None when fit was called with an array basis.                              
``icoeff``     `numpy.ndarray`_  `numpy.floating`_  Inverse-covariance diagonal, shape (nc, npoly).                                                                        
``nord``       int                                  Order of the B-spline.                                                                                                 
``npoly``      int                                  Number of polynomial basis functions in the second dimension.                                                          
``xmax``       float                                Upper normalisation bound for the polynomial basis coordinate.                                                         
``xmin``       float                                Lower normalisation bound for the polynomial basis coordinate.                                                         
=============  ================  =================  =======================================================================================================================
