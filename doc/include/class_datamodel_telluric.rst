
**Version**: 1.0.0

====================  ============================  ==========  ======================================================================================================
Attribute             Type                          Array Type  Description                                                                                           
====================  ============================  ==========  ======================================================================================================
``airmass``           float                                     Airmass of the observation                                                                            
``delta_zqso``        float                                     Allowed range for the QSO redshift about z_qso                                                        
``exptime``           float                                     Exposure time (s)                                                                                     
``func``              str                                       Polynomial function used                                                                              
``lbound_norm``       float                                     Flux normalization lower bound                                                                        
``model``             `astropy.table.table.Table`_              Table with the best-fitting model data                                                                
``npca``              int                                       Number of PCA components                                                                              
``pca_file``          str                                       Name of the QSO PCA file                                                                              
``polish``            bool                                      Perform a final optimization to tweak the best solution; see scipy.optimize.differential_evolution.   
``popsize``           int                                       A multiplier for setting the total population size for the differential evolution optimization.       
``recombination``     float                                     The recombination constant for the differential evolution optimization. Should be in the range [0, 1].
``std_cal``           str                                       File name (or shorthand) with the standard flux data                                                  
``std_dec``           float                                     DEC of the standard source                                                                            
``std_name``          str                                       Type of standard source                                                                               
``std_ra``            float                                     RA of the standard source                                                                             
``std_src``           str                                       Name of the standard source                                                                           
``telgrid``           str                                       File containing PCA components or grid of HITRAN atmosphere models                                    
``tell_norm_thresh``  float                                     ??                                                                                                    
``tell_npca``         int                                       Number of telluric PCA components used                                                                
``teltype``           str                                       Type of telluric model, `pca` or `grid`                                                               
``tol``               float                                     Relative tolerance for converage of the differential evolution optimization.                          
``ubound_norm``       float                                     Flux normalization upper bound                                                                        
``z_qso``             float                                     Redshift of the QSO                                                                                   
====================  ============================  ==========  ======================================================================================================
