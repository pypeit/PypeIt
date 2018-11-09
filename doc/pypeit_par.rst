.. highlight:: rest

.. _configobj: http://configobj.readthedocs.io/en/latest/

.. _pypitpar:

================
PypeIt Parameters
================

PypeIt allows you to customize its execution without having to change the
code directly.  While not ubiquitous, most optional arguments of pypit's
algorithms are contained within the :class:`pypit.par.pypitpar.PypeItPar`
superset.  See the `Current PypeItPar Parameter Hierarchy`_ below for the
current structure of a :class:`pypit.par.pypitpar.PypeItPar` instance.

Users can alter these parameters via the pypit file, see
:ref:`pypit_file`.  Only those parameters that the user wishes to be
different from the defaults need to be includes in the pypit file.

PypeIt has global defaults, defaults for each instrument served, and
defaults for each named pipeline approach (e.g., ARMS, etc), which are
merged in succession before starting the reduction of your data.  The
last parameters merged are those altered by the input pypit file.

PypeIt uses the `configobj`_ class to parse the user supplied arguments.
The syntax is important and the nesting of the parameter changes must
match the `Current PypeItPar Parameter Hierarchy`_.  Examples of `How to
change parameters using the pypit file`_ are given below.


Current PypeItPar Parameter Hierarchy
++++++++++++++++++++++++++++++++++++

`PypeItPar Keywords`_

    ``[rdx]``: `ReducePar Keywords`_

    ``[calibrations]``: `CalibrationsPar Keywords`_

        ``[[biasframe]]``: `FrameGroupPar Keywords`_

            ``[[[process]]]``: `ProcessImagesPar Keywords`_

        ``[[darkframe]]``: `FrameGroupPar Keywords`_

            ``[[[process]]]``: `ProcessImagesPar Keywords`_

        ``[[arcframe]]``: `FrameGroupPar Keywords`_

            ``[[[process]]]``: `ProcessImagesPar Keywords`_

        ``[[pixelflatframe]]``: `FrameGroupPar Keywords`_

            ``[[[process]]]``: `ProcessImagesPar Keywords`_

        ``[[pinholeframe]]``: `FrameGroupPar Keywords`_

            ``[[[process]]]``: `ProcessImagesPar Keywords`_

        ``[[traceframe]]``: `FrameGroupPar Keywords`_

            ``[[[process]]]``: `ProcessImagesPar Keywords`_

        ``[[standardframe]]``: `FrameGroupPar Keywords`_

            ``[[[process]]]``: `ProcessImagesPar Keywords`_

        ``[[flatfield]]``: `FlatFieldPar Keywords`_

        ``[[wavelengths]]``: `WavelengthSolutionPar Keywords`_

        ``[[slits]]``: `TraceSlitsPar Keywords`_

        ``[[tilts]]``: `WaveTiltsPar Keywords`_

    ``[scienceframe]``: `FrameGroupPar Keywords`_

        ``[[process]]``: `ProcessImagesPar Keywords`_

    ``[scienceimage]``: `ScienceImagePar Keywords`_

    ``[flexure]``: `FlexurePar Keywords`_

    ``[fluxcalib]``: `FluxCalibrationPar Keywords`_


----

PypeItPar Keywords
------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.PypeItPar`

================  ================================================  =======  ==============================  ======================================================================================================================================================================================================================================================================================
Key               Type                                              Options  Default                         Description                                                                                                                                                                                                                                                                           
================  ================================================  =======  ==============================  ======================================================================================================================================================================================================================================================================================
``rdx``           :class:`pypeit.par.pypeitpar.ReducePar`           ..       `ReducePar Keywords`_           PypIt reduction rules.                                                                                                                                                                                                                                                                
``calibrations``  :class:`pypeit.par.pypeitpar.CalibrationsPar`     ..       `CalibrationsPar Keywords`_     Parameters for the calibration algorithms                                                                                                                                                                                                                                             
``scienceframe``  :class:`pypeit.par.pypeitpar.FrameGroupPar`       ..       `FrameGroupPar Keywords`_       The frames and combination rules for the science observations                                                                                                                                                                                                                         
``scienceimage``  :class:`pypeit.par.pypeitpar.ScienceImagePar`     ..       `ScienceImagePar Keywords`_     Parameters determining sky-subtraction, object finding, and extraction                                                                                                                                                                                                                
``flexure``       :class:`pypeit.par.pypeitpar.FlexurePar`          ..       `FlexurePar Keywords`_          Parameters used by the flexure-correction procedure.  Flexure corrections are not performed by default.  To turn on, either set the parameters in the 'flexure' parameter group or set 'flexure = True' in the 'rdx' parameter group to use the default flexure-correction parameters.
``fluxcalib``     :class:`pypeit.par.pypeitpar.FluxCalibrationPar`  ..       `FluxCalibrationPar Keywords`_  Parameters used by the flux-calibration procedure.  Flux calibration is not performed by default.  To turn on, either set the parameters in the 'fluxcalib' parameter group or set 'fluxcalib = True' in the 'rdx' parameter group to use the default flux-calibration parameters.    
================  ================================================  =======  ==============================  ======================================================================================================================================================================================================================================================================================


----

ReducePar Keywords
------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.ReducePar`

================  ==========  ==============================================================================================================================================================================================================================================================================================================================================================  ============================================  ===========================================================================================================================================================================================================================================================================================================================================================
Key               Type        Options                                                                                                                                                                                                                                                                                                                                                         Default                                       Description                                                                                                                                                                                                                                                                                                                                                
================  ==========  ==============================================================================================================================================================================================================================================================================================================================================================  ============================================  ===========================================================================================================================================================================================================================================================================================================================================================
``spectrograph``  str         ``gemini_gnirs``, ``keck_deimos``, ``keck_lris_blue``, ``keck_lris_red``, ``keck_nires``, ``keck_nirspec``, ``shane_kast_blue``, ``shane_kast_red``, ``shane_kast_red_ret``, ``tng_dolores``, ``wht_isis_blue``, ``vlt_xshooter_uvb``, ``vlt_xshooter_vis``, ``vlt_xshooter_nir``, ``gemini_gmos_south``, ``gemini_gmos_north_e2v``, ``gemini_gmos_north_ham``  ..                                            Spectrograph that provided the data to be reduced.  Options are: gemini_gnirs, keck_deimos, keck_lris_blue, keck_lris_red, keck_nires, keck_nirspec, shane_kast_blue, shane_kast_red, shane_kast_red_ret, tng_dolores, wht_isis_blue, vlt_xshooter_uvb, vlt_xshooter_vis, vlt_xshooter_nir, gemini_gmos_south, gemini_gmos_north_e2v, gemini_gmos_north_ham
``detnum``        int, list   ..                                                                                                                                                                                                                                                                                                                                                              ..                                            Restrict reduction to a list of detector indices                                                                                                                                                                                                                                                                                                           
``sortroot``      str         ..                                                                                                                                                                                                                                                                                                                                                              ..                                            A filename given to output the details of the sorted files.  If None, the default is the root name of the pypeit file.  If off, no output is produced.                                                                                                                                                                                                     
``calwin``        int, float  ..                                                                                                                                                                                                                                                                                                                                                              0                                             The window of time in hours to search for calibration frames for a science frame                                                                                                                                                                                                                                                                           
``scidir``        str         ..                                                                                                                                                                                                                                                                                                                                                              ``Science``                                   Directory relative to calling directory to write science files.                                                                                                                                                                                                                                                                                            
``qadir``         str         ..                                                                                                                                                                                                                                                                                                                                                              ``QA``                                        Directory relative to calling directory to write quality assessment files.                                                                                                                                                                                                                                                                                 
``redux_path``    str         ..                                                                                                                                                                                                                                                                                                                                                              ``/Users/westfall/Work/packages/pypeit/doc``  Path to folder for performing reductions.                                                                                                                                                                                                                                                                                                                  
================  ==========  ==============================================================================================================================================================================================================================================================================================================================================================  ============================================  ===========================================================================================================================================================================================================================================================================================================================================================


----

CalibrationsPar Keywords
------------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.CalibrationsPar`

==================  ===================================================  ====================  =================================  ==================================================================================================================================================================================================
Key                 Type                                                 Options               Default                            Description                                                                                                                                                                                       
==================  ===================================================  ====================  =================================  ==================================================================================================================================================================================================
``caldir``          str                                                  ..                    ``MF``                             Directory relative to calling directory to write master files.                                                                                                                                    
``masters``         str                                                  ``reuse``, ``force``  ..                                 Treatment of master frames.  Use None to select the default behavior (which is?), 'reuse' to use any existing masters, and 'force' to __only__ use master frames.  Options are: None, reuse, force
``setup``           str                                                  ..                    ..                                 If masters='force', this is the setup name to be used: e.g., C_02_aa .  The detector number is ignored but the other information must match the Master Frames in the master frame folder.         
``trim``            bool                                                 ..                    True                               Trim the frame to isolate the data                                                                                                                                                                
``badpix``          bool                                                 ..                    True                               Make a bad pixel mask? Bias frames must be provided.                                                                                                                                              
``biasframe``       :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..                    `FrameGroupPar Keywords`_          The frames and combination rules for the bias correction                                                                                                                                          
``darkframe``       :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..                    `FrameGroupPar Keywords`_          The frames and combination rules for the dark-current correction                                                                                                                                  
``arcframe``        :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..                    `FrameGroupPar Keywords`_          The frames and combination rules for the wavelength calibration                                                                                                                                   
``pixelflatframe``  :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..                    `FrameGroupPar Keywords`_          The frames and combination rules for the field flattening                                                                                                                                         
``pinholeframe``    :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..                    `FrameGroupPar Keywords`_          The frames and combination rules for the pinholes                                                                                                                                                 
``traceframe``      :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..                    `FrameGroupPar Keywords`_          The frames and combination rules for images used for slit tracing                                                                                                                                 
``standardframe``   :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..                    `FrameGroupPar Keywords`_          The frames and combination rules for the spectrophotometric standard observations                                                                                                                 
``flatfield``       :class:`pypeit.par.pypeitpar.FlatFieldPar`           ..                    `FlatFieldPar Keywords`_           Parameters used to set the flat-field procedure                                                                                                                                                   
``wavelengths``     :class:`pypeit.par.pypeitpar.WavelengthSolutionPar`  ..                    `WavelengthSolutionPar Keywords`_  Parameters used to derive the wavelength solution                                                                                                                                                 
``slits``           :class:`pypeit.par.pypeitpar.TraceSlitsPar`          ..                    `TraceSlitsPar Keywords`_          Define how the slits should be traced using the trace ?PINHOLE? frames                                                                                                                            
``tilts``           :class:`pypeit.par.pypeitpar.WaveTiltsPar`           ..                    `WaveTiltsPar Keywords`_           Define how to tract the slit tilts using the trace frames                                                                                                                                         
==================  ===================================================  ====================  =================================  ==================================================================================================================================================================================================


----

FlatFieldPar Keywords
---------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.FlatFieldPar`

================  ====  ===========  =============  ==========================================================================================================
Key               Type  Options      Default        Description                                                                                               
================  ====  ===========  =============  ==========================================================================================================
``frame``         str   ..           ``pixelflat``  Frame to use for field flattening.  Options are: "pixelflat", or a specified calibration filename.        
``illumflatten``  bool  ..           True           Use the flat field to determine the illumination profile of each slit.                                    
``tweak_slits``   bool  ..           True           Use the illumination flat field to tweak the slit edges. illumflatten must be set to true for this to work
``method``        str   ``bspline``  ``bspline``    Method used to flat field the data; use None to skip flat-fielding.  Options are: None, bspline           
================  ====  ===========  =============  ==========================================================================================================


----

WavelengthSolutionPar Keywords
------------------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.WavelengthSolutionPar`

====================  ================  =================================  ================  =========================================================================================================================================================================================================================================================
Key                   Type              Options                            Default           Description                                                                                                                                                                                                                                              
====================  ================  =================================  ================  =========================================================================================================================================================================================================================================================
``reference``         str               ``arc``, ``sky``, ``pixel``        ``arc``           Perform wavelength calibration with an arc, sky frame.  Use 'pixel' for no wavelength solution.                                                                                                                                                          
``method``            str               ``simple``, ``fit``, ``arclines``  ``arclines``      Method to use to fit the individual arc lines.  'fit' is likely more accurate, but 'simple' uses a polynomial fit (to the log of a gaussian) and is fast and reliable.  'arclines' uses the arclines python package.Options are: simple, fit, arclines   
``lamps``             list              ..                                 ..                Name of one or more ions used for the wavelength calibration.  Use None for no calibration.  Options are: ArI, CdI, HgI, HeI, KrI, NeI, XeI, ZnI, ThAr                                                                                                   
``rms_threshold``     float             ..                                 0.15              Minimum RMS for keeping a slit solution                                                                                                                                                                                                                  
``nonlinear_counts``  float             ..                                 10000000000.0     Arc lines above this saturation threshold are not used in wavelength solution fits because they cannotbe accurately centroided                                                                                                                           
``match_toler``       float             ..                                 3.0               Matching tolerance when searching for new lines in iterative fitting of wavelength solution.This is the difference in pixels between the wavlength assigned to an arc line by an iteration of the wavelength solution to the wavelength in the line list.
``func``              str               ..                                 ``legendre``      Function used for wavelength solution fits                                                                                                                                                                                                               
``n_first``           int               ..                                 2                 Order of first guess fit to the wavelength solution.                                                                                                                                                                                                     
``n_final``           int               ..                                 4                 Order of final fit to the wavelength solution.                                                                                                                                                                                                           
``sigrej_first``      float             ..                                 2.0               Number of sigma for rejection for the first guess to the wavelength solution.                                                                                                                                                                            
``sigrej_final``      float             ..                                 3.0               Number of sigma for rejection for the final guess to the wavelength solution.                                                                                                                                                                            
``wv_cen``            float             ..                                 0.0               Central wavelength. Backwards compatibility with basic and semi-brute algorithms.                                                                                                                                                                        
``disp``              float             ..                                 0.0               Dispersion. Backwards compatibility with basic and semi-brute algorithms.                                                                                                                                                                                
``numsearch``         int               ..                                 20                Number of brightest arc lines to search for in preliminary identification                                                                                                                                                                                
``nfitpix``           int               ..                                 5                 Number of pixels to fit when deriving the centroid of the arc lines (an odd number is best)                                                                                                                                                              
``IDpixels``          int, float, list  ..                                 ..                One or more pixels at which to manually identify a line                                                                                                                                                                                                  
``IDwaves``           int, float, list  ..                                 ..                Wavelengths of the manually identified lines                                                                                                                                                                                                             
``medium``            str               ``vacuum``, ``air``                ``vacuum``        Medium used when wavelength calibrating the data.  Options are: vacuum, air                                                                                                                                                                              
``frame``             str               ``heliocentric``, ``barycentric``  ``heliocentric``  Frame of reference for the wavelength calibration.  Options are: heliocentric, barycentric                                                                                                                                                               
``min_nsig``          float             ..                                 10.0              Detection threshold for arc lines for "standard" lines                                                                                                                                                                                                   
``lowest_nsig``       float             ..                                 5.0               Detection threshold for arc lines for "weakest" lines                                                                                                                                                                                                    
====================  ================  =================================  ================  =========================================================================================================================================================================================================================================================


----

TraceSlitsPar Keywords
----------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.TraceSlitsPar`

==================  ==========  ===========================================  ============  ============================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                 Type        Options                                      Default       Description                                                                                                                                                                                                                                                                                                                                                                                                                                 
==================  ==========  ===========================================  ============  ============================================================================================================================================================================================================================================================================================================================================================================================================================================
``function``        str         ``polynomial``, ``legendre``, ``chebyshev``  ``legendre``  Function use to trace the slit center.  Options are: polynomial, legendre, chebyshev                                                                                                                                                                                                                                                                                                                                                        
``polyorder``       int         ..                                           3             Order of the function to use.                                                                                                                                                                                                                                                                                                                                                                                                               
``medrep``          int         ..                                           0             Number of times to median smooth a trace image prior to analysis for slit/order edges                                                                                                                                                                                                                                                                                                                                                       
``number``          int         ..                                           -1            Manually set the number of slits to identify (>=1). 'auto' or -1 will automatically identify the number of slits.                                                                                                                                                                                                                                                                                                                           
``trim``            tuple       ..                                           3, 3          How much to trim off each edge of each slit                                                                                                                                                                                                                                                                                                                                                                                                 
``maxgap``          int         ..                                           ..            Maximum number of pixels to allow for the gap between slits.  Use None if the neighbouring slits are far apart or of similar illumination.                                                                                                                                                                                                                                                                                                  
``maxshift``        int, float  ..                                           0.15          Maximum shift in trace crude                                                                                                                                                                                                                                                                                                                                                                                                                
``pad``             int         ..                                           0             Integer number of pixels to consider beyond the slit edges.                                                                                                                                                                                                                                                                                                                                                                                 
``sigdetect``       int, float  ..                                           20.0          Sigma detection threshold for edge detection                                                                                                                                                                                                                                                                                                                                                                                                
``fracignore``      float       ..                                           0.01          If a slit spans less than this fraction over the spectral size of the detector, it will be ignored (and reconstructed when/if an 'order' PCA analysis is performed).                                                                                                                                                                                                                                                                        
``min_slit_width``  float       ..                                           6.0           If a slit spans less than this number of arcseconds over the spatial direction of the detector, it will be ignored. Use this option to prevent the of alignment (box) slits from multislit reductions, which typically cannot be reduced without a significant struggle                                                                                                                                                                     
``diffpolyorder``   int         ..                                           2             Order of the 2D function used to fit the 2d solution for the spatial size of all orders.                                                                                                                                                                                                                                                                                                                                                    
``single``          list        ..                                           []            Add a single, user-defined slit based on its location on each detector.  Syntax is a list of values, 2 per detector, that define the slit according to column values.  The second value (for the right edge) must be greater than 0 to be applied.  LRISr example: setting single = -1, -1, 7, 295 means the code will skip the user-definition for the first detector but adds one for the second.  None means no user-level slits defined.
``sobel_mode``      str         ``nearest``, ``constant``                    ``nearest``   Mode for Sobel filtering.  Default is 'nearest' but the developers find 'constant' works best for DEIMOS.                                                                                                                                                                                                                                                                                                                                   
``pcatype``         str         ``pixel``, ``order``                         ``pixel``     Select to perform the PCA using the pixel position (pcatype=pixel) or by spectral order (pcatype=order).  Pixel positions can be used for multi-object spectroscopy where the gap between slits is irregular.  Order is used for echelle spectroscopy or for slits with separations that are a smooth function of the slit number.                                                                                                          
``pcapar``          list        ..                                           3, 2, 1, 0    Order of the polynomials to be used to fit the principle components.  The list length must be equal to or less than polyorder+1. TODO: Provide more explanation                                                                                                                                                                                                                                                                             
``pcaextrap``       list        ..                                           0, 0          The number of extra orders to predict in the negative (first number) and positive (second number) direction.  Must be two numbers in the list and they must be integers.                                                                                                                                                                                                                                                                    
==================  ==========  ===========================================  ============  ============================================================================================================================================================================================================================================================================================================================================================================================================================================


----

WaveTiltsPar Keywords
---------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.WaveTiltsPar`

===============  =========================  =============================================================  ============  ===================================================================================================================================
Key              Type                       Options                                                        Default       Description                                                                                                                        
===============  =========================  =============================================================  ============  ===================================================================================================================================
``idsonly``      bool                       ..                                                             False         Only use the arc lines that have an identified wavelength to trace tilts                                                           
``tracethresh``  int, float, list, ndarray  ..                                                             20.0          Significance threshold for arcs to be used in tracing wavelength tilts.                                                            
``order``        int                        ..                                                             2             Order of the polynomial function to be used for the tilt of an individual arc line.  Must be 1 for echelle data (Echelle pipeline).
``function``     str                        ..                                                             ``legendre``  Type of function for arc line fits                                                                                                 
``yorder``       int                        ..                                                             4             Order of the polynomial function to be used to fit the tilts along the y direction.                                                
``func2D``       str                        ..                                                             ``legendre``  Type of function for 2D fit                                                                                                        
``method``       str                        ``pca``, ``spca``, ``spline``, ``interp``, ``perp``, ``zero``  ``spca``      Method used to trace the tilt of the slit along an order.  Options are: pca, spca, spline, interp, perp, zero                      
``params``       int, list                  ..                                                             1, 1, 0       Parameters to use for the provided method.  TODO: Need more explanation                                                            
===============  =========================  =============================================================  ============  ===================================================================================================================================


----

FrameGroupPar Keywords
----------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.FrameGroupPar`

=============  ==============================================  ======================================================================================================  ============================  ===============================================================================================================================================================================================================================================================
Key            Type                                            Options                                                                                                 Default                       Description                                                                                                                                                                                                                                                    
=============  ==============================================  ======================================================================================================  ============================  ===============================================================================================================================================================================================================================================================
``frametype``  str                                             ``bias``, ``dark``, ``pixelflat``, ``arc``, ``pinhole``, ``trace``, ``standard``, ``science``, ``all``  ``bias``                      Frame type.  Options are: bias, dark, pixelflat, arc, pinhole, trace, standard, science, all                                                                                                                                                                   
``useframe``   str                                             ..                                                                                                      ..                            A master calibrations file to use if it exists.                                                                                                                                                                                                                
``number``     int                                             ..                                                                                                      0                             Used in matching calibration frames to science frames.  This sets the number of frames to use of this type                                                                                                                                                     
``exprng``     list                                            ..                                                                                                      None, None                    Used in identifying frames of this type.  This sets the minimum and maximum allowed exposure times.  There must be two items in the list.  Use None to indicate no limit; i.e., to select exposures with any time greater than 30 sec, use exprng = [30, None].
``process``    :class:`pypeit.par.pypeitpar.ProcessImagesPar`  ..                                                                                                      `ProcessImagesPar Keywords`_  Parameters used for basic image processing                                                                                                                                                                                                                     
=============  ==============================================  ======================================================================================================  ============================  ===============================================================================================================================================================================================================================================================


----

ProcessImagesPar Keywords
-------------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.ProcessImagesPar`

================  ==========  =====================================================================  ==============  ===========================================================================================================================================================================================================================================
Key               Type        Options                                                                Default         Description                                                                                                                                                                                                                                
================  ==========  =====================================================================  ==============  ===========================================================================================================================================================================================================================================
``overscan``      str         ``polynomial``, ``savgol``, ``median``                                 ``savgol``      Method used to fit the overscan.  Options are: polynomial, savgol, median                                                                                                                                                                  
``overscan_par``  int, list   ..                                                                     5, 65           Parameters for the overscan subtraction.  For 'polynomial', set overcan_par = order, number of pixels, number of repeats ; for 'savgol', set overscan_par = order, window size ; for 'median', set overscan_par = None or omit the keyword.
``match``         int, float  ..                                                                     -1              (Deprecate?) Match frames with pixel counts that are within N-sigma of one another, where match=N below.  If N < 0, nothing is matched.                                                                                                    
``combine``       str         ``mean``, ``median``, ``weightmean``                                   ``weightmean``  Method used to combine frames.  Options are: mean, median, weightmean                                                                                                                                                                      
``satpix``        str         ``reject``, ``force``, ``nothing``                                     ``reject``      Handling of saturated pixels.  Options are: reject, force, nothing                                                                                                                                                                         
``sigrej``        int, float  ..                                                                     20.0            Sigma level to reject cosmic rays (<= 0.0 means no CR removal)                                                                                                                                                                             
``n_lohi``        list        ..                                                                     0, 0            Number of pixels to reject at the lowest and highest ends of the distribution; i.e., n_lohi = low, high.  Use None for no limit.                                                                                                           
``sig_lohi``      list        ..                                                                     3.0, 3.0        Sigma-clipping level at the low and high ends of the distribution; i.e., sig_lohi = low, high.  Use None for no limit.                                                                                                                     
``replace``       str         ``min``, ``max``, ``mean``, ``median``, ``weightmean``, ``maxnonsat``  ``maxnonsat``   If all pixels are rejected, replace them using this method.  Options are: min, max, mean, median, weightmean, maxnonsat                                                                                                                    
``lamaxiter``     int         ..                                                                     1               Maximum number of iterations for LA cosmics routine.                                                                                                                                                                                       
``grow``          int, float  ..                                                                     1.5             Factor by which to expand regions with cosmic rays detected by the LA cosmics routine.                                                                                                                                                     
``rmcompact``     bool        ..                                                                     True            Remove compact detections in LA cosmics routine                                                                                                                                                                                            
``sigclip``       int, float  ..                                                                     5.0             Sigma level for rejection in LA cosmics routine                                                                                                                                                                                            
``sigfrac``       int, float  ..                                                                     0.3             Fraction for the lower clipping threshold in LA cosmics routine.                                                                                                                                                                           
``objlim``        int, float  ..                                                                     5.0             Object detection limit in LA cosmics routine                                                                                                                                                                                               
================  ==========  =====================================================================  ==============  ===========================================================================================================================================================================================================================================


----

ScienceImagePar Keywords
------------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.ScienceImagePar`

===================  ==========  =======  =======  ================================================================================
Key                  Type        Options  Default  Description                                                                     
===================  ==========  =======  =======  ================================================================================
``bspline_spacing``  int, float  ..       0.6      Break-point spacing for the bspline fit                                         
``maxnumber``        int         ..       ..       Maximum number of objects to extract in a science frame.  Use None for no limit.
``manual``           list        ..       ..       List of manual extraction parameter sets                                        
``nodding``          bool        ..       False    Use the nodded frames to perform the sky subtraction                            
===================  ==========  =======  =======  ================================================================================


----

FlexurePar Keywords
-------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.FlexurePar`

============  ==========  =================================  ==========  ====================================================================================================================================================================================================
Key           Type        Options                            Default     Description                                                                                                                                                                                         
============  ==========  =================================  ==========  ====================================================================================================================================================================================================
``method``    str         ``boxcar``, ``slitcen``, ``skip``  ``boxcar``  Method used to correct for flexure. Use None for no correction.  If slitcen is used, the flexure correction is performed before the extraction of objects.  Options are: None, boxcar, slitcen, skip
``maxshift``  int, float  ..                                 20          Maximum allowed flexure shift in pixels.                                                                                                                                                            
``spectrum``  str         ..                                 ..          Archive sky spectrum to be used for the flexure correction.                                                                                                                                         
============  ==========  =================================  ==========  ====================================================================================================================================================================================================


----

FluxCalibrationPar Keywords
---------------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.FluxCalibrationPar`

=============  ====  =======  =======  =================================================================================================================================================================
Key            Type  Options  Default  Description                                                                                                                                                      
=============  ====  =======  =======  =================================================================================================================================================================
``nonlinear``  bool  ..       False    Perform a non-linear correction.  Requires a series of pixelflats of the same lamp and setup and with a variety of exposure times and count rates in every pixel.
``sensfunc``   str   ..       ..       YAML file with an existing calibration function                                                                                                                  
=============  ====  =======  =======  =================================================================================================================================================================


