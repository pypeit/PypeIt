
**Version**: 1.2.0

====================  ============================  =================  ===============================================================================================================================================================
Attribute             Type                          Array Type         Description                                                                                                                                                    
====================  ============================  =================  ===============================================================================================================================================================
``PYP_SPEC``          str                                              PypeIt spectrograph name                                                                                                                                       
``bpmtilts``          `numpy.ndarray`_              `numpy.integer`_   Bad pixel mask for tilt solutions. Keys are taken from SlitTraceSetBitmask                                                                                     
``coeffs``            `numpy.ndarray`_              `numpy.floating`_  2D coefficents for the fit on the initial slits.  One set per slit/order (3D array).                                                                           
``func2d``            str                                              Function used for the 2D fit                                                                                                                                   
``nslit``             int                                              Total number of slits.  This can include masked slits                                                                                                          
``slits_filename``    str                                              Path to SlitTraceSet file. This helps to find the Slits calibration file when running pypeit_chk_tilts()                                                       
``spat_flexure``      float                                            Flexure shift from the input TiltImage                                                                                                                         
``spat_id``           `numpy.ndarray`_              `numpy.integer`_   Slit spat_id                                                                                                                                                   
``spat_order``        `numpy.ndarray`_              `numpy.integer`_   Order for spatial fit (nslit)                                                                                                                                  
``spec_order``        `numpy.ndarray`_              `numpy.integer`_   Order for spectral fit (nslit)                                                                                                                                 
``tilt_traces``       `astropy.table.table.Table`_                     Table with the positions of the traced and fitted tilts for all the slits. see :func:`~pypeit.wavetilts.BuildWaveTilts.make_tbl_tilt_traces` for more details. 
``tiltimg_filename``  str                                              Path to Tiltimg file. This helps to find Tiltimg file when running pypeit_chk_tilts()                                                                          
====================  ============================  =================  ===============================================================================================================================================================
