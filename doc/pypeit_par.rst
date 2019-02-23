.. highlight:: rest

.. _configobj: http://configobj.readthedocs.io/en/latest/

.. _pypeitpar:

================
PypeIt Parameters
================

PypeIt allows you to customize its execution without having to change the
code directly.

Although not ubiquitous, most optional arguments of PypeIt's
algorithms are contained within the :class:`pypeit.par.pypeitpar.PypeItPar`
superset.  See the `Current PypeItPar Parameter Hierarchy`_ below for the
current structure of a :class:`pypeit.par.pypeitpar.PypeItPar` instance.

More importantly, each instrument served provides its own default values
for :class:`pypeit.par.pypeitpar.PypeItPar` as defined by its
`default_pypeit_par` method; e.g.,
:func:`pypeit.spectrographs.shane_kast.ShaneKastSpectrograph.default_pypeit_par`.
Users can alter these parameters via the PypeIt file, see
:ref:`pypeit_file`.  Only those parameters that the user wishes to be
different from the default *used for their specified instrument* need to
be includes in the PypeIt file.

PypeIt uses the `configobj`_ class to parse the user supplied arguments.
The syntax is important and the nesting of the parameter changes must
match the `Current PypeItPar Parameter Hierarchy`_.  Examples of `How to
change parameters using the PypeIt file`_ are given below.


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

======================  ==========  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================  ========================================  ===========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                     Type        Options                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Default                                   Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
======================  ==========  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================  ========================================  ===========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``spectrograph``        str         ``gemini_gnirs``, ``keck_deimos``, ``keck_lris_blue``, ``keck_lris_red``, ``keck_nires``, ``keck_hires_red``, ``keck_hires_blue``, ``mmt_binospec``, ``keck_nirspec_low``, ``shane_kast_blue``, ``shane_kast_red``, ``shane_kast_red_ret``, ``tng_dolores``, ``wht_isis_blue``, ``vlt_xshooter_uvb``, ``vlt_xshooter_vis``, ``magellan_fire``, ``magellan_mage``, ``vlt_xshooter_nir``, ``gemini_gmos_south_ham``, ``gemini_gmos_north_e2v``, ``gemini_gmos_north_ham``, ``lbt_mods1r``, ``lbt_mods1b``, ``lbt_mods2r``, ``lbt_mods2b``, ``vlt_fors2``  ..                                        Spectrograph that provided the data to be reduced.  Options are: gemini_gnirs, keck_deimos, keck_lris_blue, keck_lris_red, keck_nires, keck_hires_red, keck_hires_blue, mmt_binospec, keck_nirspec_low, shane_kast_blue, shane_kast_red, shane_kast_red_ret, tng_dolores, wht_isis_blue, vlt_xshooter_uvb, vlt_xshooter_vis, magellan_fire, magellan_mage, vlt_xshooter_nir, gemini_gmos_south_ham, gemini_gmos_north_e2v, gemini_gmos_north_ham, lbt_mods1r, lbt_mods1b, lbt_mods2r, lbt_mods2b, vlt_fors2
``detnum``              int, list   ..                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ..                                        Restrict reduction to a list of detector indices                                                                                                                                                                                                                                                                                                                                                                                                                                                           
``sortroot``            str         ..                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ..                                        A filename given to output the details of the sorted files.  If None, the default is the root name of the pypeit file.  If off, no output is produced.                                                                                                                                                                                                                                                                                                                                                     
``calwin``              int, float  ..                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0                                         The window of time in hours to search for calibration frames for a science frame                                                                                                                                                                                                                                                                                                                                                                                                                           
``scidir``              str         ..                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ``Science``                               Directory relative to calling directory to write science files.                                                                                                                                                                                                                                                                                                                                                                                                                                            
``qadir``               str         ..                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ``QA``                                    Directory relative to calling directory to write quality assessment files.                                                                                                                                                                                                                                                                                                                                                                                                                                 
``redux_path``          str         ..                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ``/home/xavier/local/Python/PypeIt/doc``  Path to folder for performing reductions.                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``ignore_bad_headers``  bool        ..                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      False                                     Ignore bad headers (NOT recommended unless you know it is safe).                                                                                                                                                                                                                                                                                                                                                                                                                                           
======================  ==========  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================  ========================================  ===========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

CalibrationsPar Keywords
------------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.CalibrationsPar`

==================  ===================================================  =======  =================================  =========================================================================================================================================================================================
Key                 Type                                                 Options  Default                            Description                                                                                                                                                                              
==================  ===================================================  =======  =================================  =========================================================================================================================================================================================
``caldir``          str                                                  ..       ``MF``                             Directory relative to calling directory to write master files.                                                                                                                           
``reuse_masters``   bool                                                 False    ..                                 If True PypeIt will reuse existing master frames rather than recreate them. If False, it will  recreate the master frames.                                                               
``setup``           str                                                  ..       ..                                 If masters='force', this is the setup name to be used: e.g., C_02_aa .  The detector number is ignored but the other information must match the Master Frames in the master frame folder.
``trim``            bool                                                 ..       True                               Trim the frame to isolate the data                                                                                                                                                       
``badpix``          bool                                                 ..       True                               Make a bad pixel mask? Bias frames must be provided.                                                                                                                                     
``biasframe``       :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the bias correction                                                                                                                                 
``darkframe``       :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the dark-current correction                                                                                                                         
``arcframe``        :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the wavelength calibration                                                                                                                          
``pixelflatframe``  :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the field flattening                                                                                                                                
``pinholeframe``    :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the pinholes                                                                                                                                        
``traceframe``      :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for images used for slit tracing                                                                                                                        
``standardframe``   :class:`pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the spectrophotometric standard observations                                                                                                        
``flatfield``       :class:`pypeit.par.pypeitpar.FlatFieldPar`           ..       `FlatFieldPar Keywords`_           Parameters used to set the flat-field procedure                                                                                                                                          
``wavelengths``     :class:`pypeit.par.pypeitpar.WavelengthSolutionPar`  ..       `WavelengthSolutionPar Keywords`_  Parameters used to derive the wavelength solution                                                                                                                                        
``slits``           :class:`pypeit.par.pypeitpar.TraceSlitsPar`          ..       `TraceSlitsPar Keywords`_          Define how the slits should be traced using the trace ?PINHOLE? frames                                                                                                                   
``tilts``           :class:`pypeit.par.pypeitpar.WaveTiltsPar`           ..       `WaveTiltsPar Keywords`_           Define how to tract the slit tilts using the trace frames                                                                                                                                
==================  ===================================================  =======  =================================  =========================================================================================================================================================================================


----

FlatFieldPar Keywords
---------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.FlatFieldPar`

=======================  ==========  =====================  =============  =================================================================================================================================================================================================================================================
Key                      Type        Options                Default        Description                                                                                                                                                                                                                                      
=======================  ==========  =====================  =============  =================================================================================================================================================================================================================================================
``method``               str         ``bspline``, ``skip``  ``bspline``    Method used to flat field the data; use skip to skip flat-fielding.  Options are: None, bspline, skip                                                                                                                                            
``frame``                str         ..                     ``pixelflat``  Frame to use for field flattening.  Options are: "pixelflat", or a specified calibration filename.                                                                                                                                               
``illumflatten``         bool        ..                     True           Use the flat field to determine the illumination profile of each slit.                                                                                                                                                                           
``spec_samp_fine``       int, float  ..                     1.2            bspline break point spacing in units of pixels for spectral fit to flat field blaze function.                                                                                                                                                    
``spec_samp_coarse``     int, float  ..                     50.0           bspline break point spacing in units of pixels for 2-d bspline-polynomial fit to flat field image residuals. This should be a large number unless you are trying to fit a sky flat with lots of narrow spectral features.                        
``spat_samp``            int, float  ..                     5.0            Spatial sampling for slit illumination function. This is the width of the median filter in pixels used to determine the slit illumination function, and thus sets the minimum scale on which the illumination function will have features.       
``tweak_slits``          bool        ..                     True           Use the illumination flat field to tweak the slit edges. This will work even if illumflatten is set to False                                                                                                                                     
``tweak_slits_thresh``   float       ..                     0.93           If tweak_slits is True, this sets the illumination function threshold used to tweak the slit boundaries based on the illumination flat. It should be a number less than 1.0                                                                      
``tweak_slits_maxfrac``  float       ..                     0.1            If tweak_slit is True, this sets the maximum fractional amount (of a slits width) allowed for trimming each (i.e. left and right) slit boundary, i.e. the default is 10% which means slits would shrink or grow by at most 20% (10% on each side)
=======================  ==========  =====================  =============  =================================================================================================================================================================================================================================================


----

WavelengthSolutionPar Keywords
------------------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.WavelengthSolutionPar`

====================  =========================  ========================================================================================  ================  ====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                   Type                       Options                                                                                   Default           Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
====================  =========================  ========================================================================================  ================  ====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``reference``         str                        ``arc``, ``sky``, ``pixel``                                                               ``arc``           Perform wavelength calibration with an arc, sky frame.  Use 'pixel' for no wavelength solution.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``method``            str                        ``simple``, ``semi-brute``, ``basic``, ``holy-grail``, ``reidentify``, ``full_template``  ``holy-grail``    Method to use to fit the individual arc lines. Most of these methods are now deprecated as they fail most of the time without significant parameter tweaking. 'holy-grail' attempts to get a first guess at line IDs by looking for patterns in the line locations. It is fully automated and works really well excpet for when it does not'reidentify' is now the preferred method, however it requires that an archive of wavelength solution has been constructed for your instrument/grating combination                           Options are: simple, semi-brute, basic, holy-grail, reidentify, full_template
``echelle``           bool                       ..                                                                                        False             Is this an echelle spectrograph? If yes an additional 2-d fit wavelength fit will be performed as a function of spectral pixel and order number to improve the wavelength solution                                                                                                                                                                                                                                                                                                                                                                                                                                  
``ech_fix_format``    bool                       ..                                                                                        True              Is this a fixed format echelle like ESI, X-SHOOTER, or NIRES. If so reidentification will assume that each order in the data is aligned with a single order in the reid arxiv                                                                                                                                                                                                                                                                                                                                                                                                                                       
``ech_nspec_coeff``   int                        ..                                                                                        4                 For echelle spectrographs, order of the final 2d fit to the spectral dimension. You should choose this to be the n_final of the fits to the individual orders.                                                                                                                                                                                                                                                                                                                                                                                                                                                      
``ech_norder_coeff``  int                        ..                                                                                        4                 For echelle spectrographs, order of the final 2d fit to the order dimension.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``ech_sigrej``        int, float                 ..                                                                                        2.0               For echelle spectrographs sigma clipping rejection threshold in 2d fit to spectral and order dimensions                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``lamps``             list                       ..                                                                                        ..                Name of one or more ions used for the wavelength calibration.  Use None for no calibration.  Options are: ArI, CdI, HgI, HeI, KrI, NeI, XeI, ZnI, ThAr                                                                                                                                                                                                                                                                                                                                                                                                                                                              
``nonlinear_counts``  float                      ..                                                                                        10000000000.0     Arc lines above this saturation threshold are not used in wavelength solution fits because they cannotbe accurately centroided                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
``sigdetect``         int, float, list, ndarray  ..                                                                                        5.0               Detection threshold for arc lines. This can be a single number or a list/array providing the value for each slit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``fwhm``              int, float                 ..                                                                                        4.0               Spectral sampling of the arc lines. This is the FWHM of an arcline in *unbinned* pixels.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
``reid_arxiv``        str                        ..                                                                                        ..                Name of the archival wavelength solution file that will be used for the wavelength reidentification if the wavelength solution method = reidentify                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``nreid_min``         int                        ..                                                                                        1                 Minimum number of times that a given candidate reidentified line must be properly matched with a line in the arxiv to be considered a good reidentification. If there is a lot of duplication in the arxiv of the spectra in question (i.e. multislit) set this to a number like 1-4. For echelle this depends on the number of solutions in the arxiv. For fixed format echelle (ESI, X-SHOOTER, NIRES) set this 1. For an echelle with a tiltable grating, it will depend on the number of solutions in the arxiv.                                                                                                
``cc_thresh``         float, list, ndarray       ..                                                                                        0.7               Threshold for the *global* cross-correlation coefficient between an input spectrum and member of the archive required to attempt reidentification. Spectra from the archive with a lower cross-correlation are not used for reidentification. This can be a single number or a list/array providing the value for each slit                                                                                                                                                                                                                                                                                         
``cc_local_thresh``   float                      ..                                                                                        0.7               Threshold for the *local* cross-correlation coefficient, evaluated at each reidentified line,  between an input spectrum and the shifted and stretched archive spectrum above which a line must be to be considered a good line for reidentification. The local cross-correlation is evaluated at each candidate reidentified line (using a window of nlocal_cc), and is then used to score the the reidentified lines to arrive at the final set of good reidentifications                                                                                                                                         
``nlocal_cc``         int                        ..                                                                                        11                Size of pixel window used for local cross-correlation computation for each arc line. If not an odd number one will be added to it to make it odd.                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``rms_threshold``     float, list, ndarray       ..                                                                                        0.15              Minimum RMS for keeping a slit/order solution. This can be a single number or a list/array providing the value for each slit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``match_toler``       float                      ..                                                                                        2.0               Matching tolerance in pixels when searching for new lines. This is the difference in pixels between the wavlength assigned to an arc line by an iteration of the wavelength solution to the wavelength in the line list. This parameter is also used as the matching tolerance in pixels for a line reidentification. A good line match must match within this tolerance to the shifted and stretched archive spectrum, and the archive wavelength solution at this match must be within match_toler dispersion elements from the line in line list.                                                                
``func``              str                        ..                                                                                        ``legendre``      Function used for wavelength solution fits                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
``n_first``           int                        ..                                                                                        2                 Order of first guess fit to the wavelength solution.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
``n_final``           int, float, list, ndarray  ..                                                                                        4                 Order of final fit to the wavelength solution. This can be a single number or a list/array providing the value for each slit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``sigrej_first``      float                      ..                                                                                        2.0               Number of sigma for rejection for the first guess to the wavelength solution.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
``sigrej_final``      float                      ..                                                                                        3.0               Number of sigma for rejection for the final guess to the wavelength solution.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
``wv_cen``            float                      ..                                                                                        0.0               Central wavelength. Backwards compatibility with basic and semi-brute algorithms.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``disp``              float                      ..                                                                                        0.0               Dispersion. Backwards compatibility with basic and semi-brute algorithms.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
``numsearch``         int                        ..                                                                                        20                Number of brightest arc lines to search for in preliminary identification                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
``nfitpix``           int                        ..                                                                                        5                 Number of pixels to fit when deriving the centroid of the arc lines (an odd number is best)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
``IDpixels``          int, float, list           ..                                                                                        ..                One or more pixels at which to manually identify a line                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``IDwaves``           int, float, list           ..                                                                                        ..                Wavelengths of the manually identified lines                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``medium``            str                        ``vacuum``, ``air``                                                                       ``vacuum``        Medium used when wavelength calibrating the data.  Options are: vacuum, air                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
``frame``             str                        ``heliocentric``, ``barycentric``                                                         ``heliocentric``  Frame of reference for the wavelength calibration.  Options are: heliocentric, barycentric                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
``nsnippet``          int                        ..                                                                                        2                 Number of spectra to chop the arc spectrum into when using the full_template method                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
====================  =========================  ========================================================================================  ================  ====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

TraceSlitsPar Keywords
----------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.TraceSlitsPar`

====================  ==========  ===========================================  ============  ============================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                   Type        Options                                      Default       Description                                                                                                                                                                                                                                                                                                                                                                                                                                 
====================  ==========  ===========================================  ============  ============================================================================================================================================================================================================================================================================================================================================================================================================================================
``function``          str         ``polynomial``, ``legendre``, ``chebyshev``  ``legendre``  Function use to trace the slit center.  Options are: polynomial, legendre, chebyshev                                                                                                                                                                                                                                                                                                                                                        
``polyorder``         int         ..                                           3             Order of the function to use.                                                                                                                                                                                                                                                                                                                                                                                                               
``medrep``            int         ..                                           0             Number of times to median smooth a trace image prior to analysis for slit/order edges                                                                                                                                                                                                                                                                                                                                                       
``number``            int         ..                                           -1            Manually set the number of slits to identify (>=1). 'auto' or -1 will automatically identify the number of slits.                                                                                                                                                                                                                                                                                                                           
``trim``              tuple       ..                                           0, 0          How much to trim off each edge of each slit.  Each number should be 0 or positive                                                                                                                                                                                                                                                                                                                                                           
``maxgap``            int         ..                                           ..            Maximum number of pixels to allow for the gap between slits.  Use None if the neighbouring slits are far apart or of similar illumination.                                                                                                                                                                                                                                                                                                  
``maxshift``          int, float  ..                                           0.15          Maximum shift in trace crude. Use a larger number for more curved slits/orders.                                                                                                                                                                                                                                                                                                                                                             
``pad``               int         ..                                           0             Integer number of pixels to consider beyond the slit edges.                                                                                                                                                                                                                                                                                                                                                                                 
``sigdetect``         int, float  ..                                           20.0          Sigma detection threshold for edge detection                                                                                                                                                                                                                                                                                                                                                                                                
``min_slit_width``    float       ..                                           6.0           If a slit spans less than this number of arcseconds over the spatial direction of the detector, it will be ignored. Use this option to prevent the of alignment (box) slits from multislit reductions, which typically cannot be reduced without a significant struggle                                                                                                                                                                     
``add_slits``         str, list   ..                                           []            Add one or more user-defined slits.  This is a list of lists, with each sub-list having syntax (all integers):  det:spat0:spat1:spec  where det=detector, spat=spatial pixel, spec=spectral pixel, For example,  2:2121:2322:2000,3:1201:1500:2000                                                                                                                                                                                          
``rm_slits``          str, list   ..                                           []            Remove one or more user-specified slits.  This is a list of lists, with each sub-list having syntax (all integers):  det:spat:spec where det=detector, spat=spatial pixel, spec=spectral pixel, for example,  2:2121:2000,3:1500:2000the slit tracing code will remove the slits on detector 2 that contain pixel (spat,spec)=(2121,2000)                                                                                                   
``diffpolyorder``     int         ..                                           2             Order of the 2D function used to fit the 2d solution for the spatial size of all orders.                                                                                                                                                                                                                                                                                                                                                    
``single``            list        ..                                           []            Add a single, user-defined slit based on its location on each detector.  Syntax is a list of values, 2 per detector, that define the slit according to column values.  The second value (for the right edge) must be greater than 0 to be applied.  LRISr example: setting single = -1, -1, 7, 295 means the code will skip the user-definition for the first detector but adds one for the second.  None means no user-level slits defined.
``sobel_mode``        str         ``nearest``, ``constant``                    ``nearest``   Mode for Sobel filtering.  Default is 'nearest' but the developers find 'constant' works best for DEIMOS.                                                                                                                                                                                                                                                                                                                                   
``pcatype``           str         ``pixel``, ``order``                         ``pixel``     Select to perform the PCA using the pixel position (pcatype=pixel) or by spectral order (pcatype=order).  Pixel positions can be used for multi-object spectroscopy where the gap between slits is irregular.  Order is used for echelle spectroscopy or for slits with separations that are a smooth function of the slit number.                                                                                                          
``pcapar``            list        ..                                           3, 2, 1, 0    Order of the polynomials to be used to fit the principle components.  The list length must be equal to or less than polyorder+1. TODO: Provide more explanation                                                                                                                                                                                                                                                                             
``pcaextrap``         list        ..                                           0, 0          The number of extra orders to predict in the negative (first number) and positive (second number) direction.  Must be two numbers in the list and they must be integers.                                                                                                                                                                                                                                                                    
``smash_range``       list        ..                                           0.0, 1.0      Range of the slit in the spectral direction (in fractional units) to smash when searching for slit edges. If the spectrum covers only a portion of the image, use that range.                                                                                                                                                                                                                                                               
``mask_frac_thresh``  float       ..                                           0.6           Minimum fraction of the slit edge that was *not* masked to use in initial PCA.                                                                                                                                                                                                                                                                                                                                                              
====================  ==========  ===========================================  ============  ============================================================================================================================================================================================================================================================================================================================================================================================================================================


----

WaveTiltsPar Keywords
---------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.WaveTiltsPar`

===================  =========================  =======  ==============  =========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                  Type                       Options  Default         Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
===================  =========================  =======  ==============  =========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``idsonly``          bool                       ..       False           Only use the arc lines that have an identified wavelength to trace tilts                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
``tracethresh``      int, float, list, ndarray  ..       20.0            Significance threshold for arcs to be used in tracing wavelength tilts. This can be a single number or a list/array providing the value for each slit                                                                                                                                                                                                                                                                                                                                                                                                                    
``sig_neigh``        int, float                 ..       10.0            Significance threshold for arcs to be used in line identification for the purpose of identifying neighboring lines.The tracethresh parameter above determines the significance threshold of lines that will be traced, but these lines must be at least nfwhm_neigh fwhm away from neighboring lines. This parameter determines the significance above which a line must be to be considered a possible colliding neighbor. A low value of sig_neigh will result in an overall larger number of lines, which will result in more lines above tracethresh getting rejected
``nfwhm_neigh``      int, float                 ..       3.0             Required separation between neighboring arc lines for them to be considered for tilt tracing in units of the the spectral fwhm (see wavelength parset where fwhm is defined)                                                                                                                                                                                                                                                                                                                                                                                             
``maxdev_tracefit``  int, float                 ..       0.2             Maximum absolute deviation (in units of fwhm) for the legendre polynomial fits to individual arc line tilt fits during iterative trace fitting (flux weighted, then gaussian weighted)                                                                                                                                                                                                                                                                                                                                                                                   
``sigrej_trace``     int, float                 ..       3.0             Outlier rejection significance to determine which traced arc lines should be included in the global fit                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``spat_order``       int, float, list, ndarray  ..       3               Order of the legendre polynomial to be fit to the the tilt of an arc line. This parameter determinesboth the orer of the *individual* arc line tilts, as well as the order of the spatial direction of the2d legendre polynomial (spatial, spectral) that is fit to obtain a global solution for the tilts across theslit/order. This can be a single number or a list/array providing the value for each slit                                                                                                                                                           
``spec_order``       int, float, list, ndarray  ..       4               Order of the spectral direction of the 2d legendre polynomial (spatial, spectral) that is fit to obtain a global solution for the tilts across the slit/order. This can be a single number or a list/array providing the value for each slit                                                                                                                                                                                                                                                                                                                             
``func2d``           str                        ..       ``legendre2d``  Type of function for 2D fit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
``maxdev2d``         int, float                 ..       0.25            Maximum absolute deviation (in units of fwhm) rejection threshold used to determines which pixels in global 2d fits to arc line tilts are rejected because they deviate from the model by more than this value                                                                                                                                                                                                                                                                                                                                                           
``sigrej2d``         int, float                 ..       3.0             Outlier rejection significance determining which pixels on a fit to an arc line tilt are rejected by the global 2D fit                                                                                                                                                                                                                                                                                                                                                                                                                                                   
===================  =========================  =======  ==============  =========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

FrameGroupPar Keywords
----------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.FrameGroupPar`

=============  ==============================================  ======================================================================================================  ============================  ===============================================================================================================================================================================================================================================================
Key            Type                                            Options                                                                                                 Default                       Description                                                                                                                                                                                                                                                    
=============  ==============================================  ======================================================================================================  ============================  ===============================================================================================================================================================================================================================================================
``frametype``  str                                             ``bias``, ``dark``, ``pixelflat``, ``arc``, ``pinhole``, ``trace``, ``standard``, ``science``, ``all``  ``science``                   Frame type.  Options are: bias, dark, pixelflat, arc, pinhole, trace, standard, science, all                                                                                                                                                                   
``useframe``   str                                             ..                                                                                                      ``science``                   A master calibrations file to use if it exists.                                                                                                                                                                                                                
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
``sigclip``       int, float  ..                                                                     4.5             Sigma level for rejection in LA cosmics routine                                                                                                                                                                                            
``sigfrac``       int, float  ..                                                                     0.3             Fraction for the lower clipping threshold in LA cosmics routine.                                                                                                                                                                           
``objlim``        int, float  ..                                                                     3.0             Object detection limit in LA cosmics routine                                                                                                                                                                                               
================  ==========  =====================================================================  ==============  ===========================================================================================================================================================================================================================================


----

ScienceImagePar Keywords
------------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.ScienceImagePar`

===================  ==========  =======  =======  =============================================================================================================================================================================================================================================================================================
Key                  Type        Options  Default  Description                                                                                                                                                                                                                                                                                  
===================  ==========  =======  =======  =============================================================================================================================================================================================================================================================================================
``bspline_spacing``  int, float  ..       0.6      Break-point spacing for the bspline sky subtraction fits.                                                                                                                                                                                                                                    
``sig_thresh``       int, float  ..       10.0     Significance threshold for object finding.                                                                                                                                                                                                                                                   
``maxnumber``        int         ..       10       Maximum number of objects to extract in a science frame.  Use None for no limit.                                                                                                                                                                                                             
``sn_gauss``         int, float  ..       4.0      S/N threshold for performing the more sophisticated optimal extraction which performs a b-spline fit to the object profile. For S/N < sn_gauss the code will simply optimal extractwith a Gaussian with FWHM determined from the object finding.                                             
``model_full_slit``  bool        ..       False    If True local sky subtraction will be performed on the entire slit. If False, local sky subtraction will be applied to only a restricted region around each object. This should be set to True for either multislit observations using narrow slits or echelle observations with narrow slits
``no_poly``          bool        ..       False    Turn off polynomial basis (Legendre) in global sky subtraction                                                                                                                                                                                                                               
``manual``           list        ..       ..       List of manual extraction parameter sets                                                                                                                                                                                                                                                     
===================  ==========  =======  =======  =============================================================================================================================================================================================================================================================================================


----

FlexurePar Keywords
-------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.FlexurePar`

============  ==========  =================================  ==========================================================================  ======================================================================================================================================================================================================================
Key           Type        Options                            Default                                                                     Description                                                                                                                                                                                                           
============  ==========  =================================  ==========================================================================  ======================================================================================================================================================================================================================
``method``    str         ``boxcar``, ``slitcen``, ``skip``  ``skip``                                                                    Method used to correct for flexure. Use skip for no correction.  If slitcen is used, the flexure correction is performed before the extraction of objects (not recommended).  Options are: None, boxcar, slitcen, skip
``maxshift``  int, float  ..                                 20                                                                          Maximum allowed flexure shift in pixels.                                                                                                                                                                              
``spectrum``  str         ..                                 ``/home/xavier/local/Python/PypeIt/pypeit/data/sky_spec/paranal_sky.fits``  Archive sky spectrum to be used for the flexure correction.                                                                                                                                                           
============  ==========  =================================  ==========================================================================  ======================================================================================================================================================================================================================


----

FluxCalibrationPar Keywords
---------------------------

Class Instantiation: :class:`pypeit.par.pypeitpar.FluxCalibrationPar`

=================  ========  =======  =======  =========================================================================================================================================================================================================================================
Key                Type      Options  Default  Description                                                                                                                                                                                                                              
=================  ========  =======  =======  =========================================================================================================================================================================================================================================
``balm_mask_wid``  float     ..       5.0      Mask width for Balmer lines in Angstroms.                                                                                                                                                                                                
``std_file``       str       ..       ..       Standard star file to generate sensfunc                                                                                                                                                                                                  
``std_obj_id``     str, int  ..       ..       Specifies object in spec1d file to use as standard. The brightest object found is used otherwise.                                                                                                                                        
``sensfunc``       str       ..       ..       FITS file that contains or will contain the sensitivity function.                                                                                                                                                                        
``star_type``      float     ..       ..       Spectral type of the standard star (for near-IR mainly)                                                                                                                                                                                  
``star_mag``       float     ..       ..       Magnitude of the standard star (for near-IR mainly)                                                                                                                                                                                      
``multi_det``      list      ..       ..       List of detector numbers to splice together for multi-detector instruments (e.g. DEIMOS) They are assumed to be in order of increasing wavelength And that there is *no* overlap in wavelength across detectors (might be ok if there is)
``telluric``       bool      ..       False    If telluric=True the code creates a sintetic standard star spectrum using the Kurucz models, the sens func is created setting nresln=1.5 it contains the correction for telluric lines.                                                  
=================  ========  =======  =======  =========================================================================================================================================================================================================================================



Instrument-Specific Default Configuration
+++++++++++++++++++++++++++++++++++++++++

The following provides the changes to the global default parameters
provided above for each instrument.  That is, if one were to include
these in the PypeIt file, you would be reproducing the effect of the
`default_pypeit_par` method specific to each derived
:class:`pypeit.spectrographs.spectrograph.Spectrograph` class.

KECK DEIMOS
-----------
Alterations to the default parameters are::

  [rdx]
      spectrograph = keck_deimos
  [calibrations]
      [[biasframe]]
          useframe = overscan
          number = 5
          exprng = None, 2
      [[darkframe]]
          exprng = 999999, None
      [[arcframe]]
          number = 1
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          exprng = None, 30
          [[[process]]]
              combine = median
              sig_lohi = 10.0, 10.0
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 3
          exprng = None, 30
      [[standardframe]]
          number = 1
      [[wavelengths]]
          lamps = ArI, NeI, KrI, XeI
          nonlinear_counts = 56360.1
          match_toler = 2.5
          n_first = 3
      [[slits]]
          sigdetect = 50.0
  [scienceframe]
      exprng = 30, None
      [[process]]
          sigclip = 4.0
          objlim = 1.5
  [flexure]
      method = boxcar

KECK LRISb
----------
Alterations to the default parameters are::

  [rdx]
      spectrograph = keck_lris_blue
  [calibrations]
      [[biasframe]]
          number = 5
          exprng = None, 1
      [[darkframe]]
          exprng = 999999, None
      [[arcframe]]
          number = 1
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          exprng = None, 30
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 3
          exprng = None, 30
      [[standardframe]]
          number = 1
      [[wavelengths]]
          method = full_template
          lamps = NeI, ArI, CdI, KrI, XeI, ZnI, HgI
          nonlinear_counts = 56360.1
          sigdetect = 10.0
          rms_threshold = 0.2
          match_toler = 2.5
          n_first = 3
      [[slits]]
          sigdetect = 30.0
  [scienceframe]
      exprng = 29, None
  [flexure]
      method = boxcar

KECK LRISr
----------
Alterations to the default parameters are::

  [rdx]
      spectrograph = keck_lris_red
  [calibrations]
      [[biasframe]]
          number = 5
          exprng = None, 1
      [[darkframe]]
          exprng = 999999, None
      [[arcframe]]
          number = 1
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          exprng = None, 30
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 3
          exprng = None, 30
      [[standardframe]]
          number = 1
      [[wavelengths]]
          lamps = NeI, ArI, CdI, KrI, XeI, ZnI, HgI
          nonlinear_counts = 49806.6
          sigdetect = 10.0
          rms_threshold = 0.2
      [[slits]]
          sigdetect = 50.0
      [[tilts]]
          tracethresh = 25
          maxdev_tracefit = 1.0
          spat_order = 4
          spec_order = 7
          maxdev2d = 1.0
          sigrej2d = 5.0
  [scienceframe]
      exprng = 29, None
      [[process]]
          sigclip = 5.0
          objlim = 5.0
  [scienceimage]
      bspline_spacing = 0.8
  [flexure]
      method = boxcar

KECK NIRES
----------
Alterations to the default parameters are::

  [rdx]
      spectrograph = keck_nires
  [calibrations]
      [[biasframe]]
          useframe = overscan
      [[darkframe]]
          exprng = 20, None
      [[arcframe]]
          number = 1
          exprng = 20, None
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
      [[traceframe]]
          number = 5
      [[standardframe]]
          number = 1
          exprng = None, 20
      [[flatfield]]
          illumflatten = False
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_norder_coeff = 6
          ech_sigrej = 3.0
          lamps = OH_NIRES
          nonlinear_counts = 760000.0
          fwhm = 5.0
          reid_arxiv = keck_nires.json
          rms_threshold = 0.2
          n_final = 3, 4, 4, 4, 4
      [[tilts]]
          tracethresh = 10.0
  [scienceframe]
      useframe = overscan
      exprng = 20, None
      [[process]]
          satpix = nothing
          sigclip = 20.0
  [scienceimage]
      bspline_spacing = 0.8

KECK NIRSPEC
------------
Alterations to the default parameters are::

  [calibrations]
      [[biasframe]]
          exprng = None, 2
      [[darkframe]]
          exprng = None, 5
      [[arcframe]]
          number = 1
          exprng = 1, None
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          exprng = 0, None
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 5
          exprng = 0, None
      [[standardframe]]
          number = 1
          exprng = None, 5
      [[wavelengths]]
          lamps = OH_R24000
          rms_threshold = 0.2
      [[tilts]]
          tracethresh = 10.0
  [scienceframe]
      exprng = 1, None

SHANE KASTb
-----------
Alterations to the default parameters are::

  [rdx]
      spectrograph = shane_kast_blue
  [calibrations]
      [[biasframe]]
          number = 5
          exprng = None, 1
      [[darkframe]]
          exprng = 999999, None
      [[arcframe]]
          number = 1
          exprng = None, 61
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          exprng = 0, None
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 5
          exprng = 0, None
      [[standardframe]]
          number = 1
          exprng = 1, 61
      [[wavelengths]]
          method = full_template
          lamps = CdI, HgI, HeI
          nonlinear_counts = 49806.6
          rms_threshold = 0.2
          match_toler = 2.5
          n_first = 3
      [[tilts]]
          maxdev_tracefit = 0.02
          spec_order = 5
          maxdev2d = 0.02
  [scienceframe]
      exprng = 61, None
  [flexure]
      method = boxcar
      spectrum = /home/xavier/local/Python/PypeIt/pypeit/data/sky_spec/sky_kastb_600.fits

SHANE KASTr
-----------
Alterations to the default parameters are::

  [rdx]
      spectrograph = shane_kast_red
  [calibrations]
      [[biasframe]]
          number = 5
          exprng = None, 1
      [[darkframe]]
          exprng = 999999, None
      [[arcframe]]
          number = 1
          exprng = None, 61
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          exprng = 0, None
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 5
          exprng = 0, None
      [[standardframe]]
          number = 1
          exprng = 1, 61
      [[wavelengths]]
          lamps = NeI, HgI, HeI, ArI
          nonlinear_counts = 49806.6
  [scienceframe]
      exprng = 61, None
  [flexure]
      method = boxcar

SHANE KASTr
-----------
Alterations to the default parameters are::

  [rdx]
      spectrograph = shane_kast_red_ret
  [calibrations]
      [[biasframe]]
          number = 5
          exprng = None, 1
      [[darkframe]]
          exprng = 999999, None
      [[arcframe]]
          number = 1
          exprng = None, 61
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 3
          exprng = 0, None
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 3
          exprng = 0, None
      [[standardframe]]
          number = 1
          exprng = 1, 61
      [[wavelengths]]
          lamps = NeI, HgI, HeI, ArI
          nonlinear_counts = 91200.0
  [scienceframe]
      exprng = 61, None
  [flexure]
      method = boxcar

TNG DOLORES
-----------
Alterations to the default parameters are::

  [calibrations]
      [[biasframe]]
          number = 5
          exprng = None, 0.1
      [[darkframe]]
          exprng = 999999, None
      [[arcframe]]
          number = 1
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 3
      [[standardframe]]
          number = 1
  [scienceframe]
      exprng = 1, None

WHT ISISb
---------
Alterations to the default parameters are::

  [rdx]
      spectrograph = wht_isis_blue
  [calibrations]
      [[biasframe]]
          number = 5
          exprng = None, 1
      [[darkframe]]
          exprng = 999999, None
      [[arcframe]]
          number = 1
          exprng = None, 120
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          [[[process]]]
              combine = median
              sig_lohi = 10.0, 10.0
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 3
      [[standardframe]]
          number = 1
          exprng = None, 120
      [[wavelengths]]
          method = simple
  [scienceframe]
      exprng = 90, None

VLT XShooter_UVB
----------------
Alterations to the default parameters are::

  [rdx]
      spectrograph = vlt_xshooter_uvb
  [calibrations]
      [[biasframe]]
          number = 5
      [[arcframe]]
          number = 1
          [[[process]]]
              overscan = median
              sigrej = -1
      [[pixelflatframe]]
          number = 5
      [[traceframe]]
          number = 3
          [[[process]]]
              overscan = median
      [[standardframe]]
          number = 1
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_norder_coeff = 5
          ech_sigrej = 3.0
          lamps = ThAr_XSHOOTER_UVB
          nonlinear_counts = 55900.0
          reid_arxiv = vlt_xshooter_uvb1x1_iraf.json
          rms_threshold = 0.5
      [[slits]]
          polyorder = 5
          maxshift = 0.5
          sigdetect = 8.0
  [scienceframe]
      useframe = overscan

VLT XShooter_VIS
----------------
Alterations to the default parameters are::

  [rdx]
      spectrograph = vlt_xshooter_vis
  [calibrations]
      [[biasframe]]
          useframe = overscan
          number = 5
      [[arcframe]]
          useframe = overscan
          number = 1
          [[[process]]]
              overscan = median
              sigrej = -1
      [[pixelflatframe]]
          number = 5
      [[traceframe]]
          useframe = overscan
          number = 3
          [[[process]]]
              overscan = median
      [[standardframe]]
          number = 1
      [[flatfield]]
          tweak_slits_thresh = 0.9
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_sigrej = 3.0
          lamps = ThAr_XSHOOTER_VIS
          nonlinear_counts = 56360.1
          fwhm = 11.0
          reid_arxiv = vlt_xshooter_vis1x1.json
          cc_thresh = 0.5
          cc_local_thresh = 0.5
          rms_threshold = 0.5
          n_final = 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3
      [[slits]]
          polyorder = 6
          maxshift = 0.5
          sigdetect = 8.0
      [[tilts]]
          tracethresh = 15
          spec_order = 5
  [scienceframe]
      useframe = overscan
  [scienceimage]
      bspline_spacing = 0.8
      model_full_slit = True

VLT XShooter_NIR
----------------
Alterations to the default parameters are::

  [rdx]
      spectrograph = vlt_xshooter_nir
  [calibrations]
      [[biasframe]]
          useframe = none
          number = 5
      [[arcframe]]
          number = 1
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
      [[traceframe]]
          number = 3
      [[standardframe]]
          number = 1
      [[flatfield]]
          illumflatten = False
          tweak_slits_thresh = 0.9
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_nspec_coeff = 5
          ech_norder_coeff = 5
          ech_sigrej = 3.0
          lamps = OH_XSHOOTER
          nonlinear_counts = 172000.0
          sigdetect = 10.0
          fwhm = 5.0
          reid_arxiv = vlt_xshooter_nir.json
          rms_threshold = 0.25
      [[slits]]
          polyorder = 5
          maxshift = 0.5
          sigdetect = 120.0
          pcatype = order
      [[tilts]]
          tracethresh = 25.0
          maxdev_tracefit = 0.04
          maxdev2d = 0.04
  [scienceframe]
      useframe = none
      [[process]]
          satpix = nothing
          sigclip = 20.0
  [scienceimage]
      bspline_spacing = 0.8
      model_full_slit = True

GEMINI-N GNIRS
--------------
Alterations to the default parameters are::

  [rdx]
      spectrograph = gemini_gnirs
  [calibrations]
      [[biasframe]]
          useframe = overscan
      [[arcframe]]
          number = 1
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          exprng = None, 30
      [[traceframe]]
          number = 5
          exprng = None, 30
      [[standardframe]]
          number = 1
          exprng = None, 30
      [[flatfield]]
          illumflatten = False
          tweak_slits_thresh = 0.9
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_nspec_coeff = 3
          ech_norder_coeff = 5
          ech_sigrej = 3.0
          lamps = OH_GNIRS
          nonlinear_counts = 63900.0
          reid_arxiv = gemini_gnirs.json
          cc_thresh = 0.6
          rms_threshold = 1.0
          n_final = 1, 3, 3, 3, 3, 3
      [[slits]]
          polyorder = 5
          maxshift = 0.5
          sigdetect = 50.0
      [[tilts]]
          tracethresh = 5.0, 10, 10, 10, 10, 10
          sig_neigh = 5.0
          nfwhm_neigh = 2.0
  [scienceframe]
      useframe = overscan
      exprng = 30, None
  [scienceimage]
      bspline_spacing = 0.8
      sig_thresh = 5.0
      model_full_slit = True

GEMINI-S GMOS-S
---------------
Alterations to the default parameters are::

  [calibrations]
      [[biasframe]]
          number = 5
      [[arcframe]]
          number = 1
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          [[[process]]]
              combine = median
              sig_lohi = 10.0, 10.0
      [[traceframe]]
          number = 3
      [[standardframe]]
          number = 1
      [[wavelengths]]
          lamps = CuI, ArI, ArII
          rms_threshold = 0.4

GEMINI-N GMOS-N
---------------
Alterations to the default parameters are::

  [calibrations]
      [[biasframe]]
          number = 5
      [[arcframe]]
          number = 1
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          [[[process]]]
              combine = median
              sig_lohi = 10.0, 10.0
      [[traceframe]]
          number = 3
      [[standardframe]]
          number = 1
      [[wavelengths]]
          lamps = CuI, ArI, ArII
          rms_threshold = 0.4

GEMINI-N GMOS-N
---------------
Alterations to the default parameters are::

  [calibrations]
      [[biasframe]]
          number = 5
      [[arcframe]]
          number = 1
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          [[[process]]]
              combine = median
              sig_lohi = 10.0, 10.0
      [[traceframe]]
          number = 3
      [[standardframe]]
          number = 1
      [[wavelengths]]
          lamps = CuI, ArI, ArII
          rms_threshold = 0.4

MAGELLAN FIRE
-------------
Alterations to the default parameters are::

  [rdx]
      spectrograph = magellan_fire
  [calibrations]
      [[biasframe]]
          useframe = overscan
      [[darkframe]]
          exprng = 20, None
      [[arcframe]]
          number = 1
          exprng = 20, None
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
      [[traceframe]]
          number = 5
      [[standardframe]]
          number = 1
          exprng = None, 60
      [[wavelengths]]
          echelle = True
          ech_sigrej = 3.0
          lamps = OH_XSHOOTER
          nonlinear_counts = 20000.0
          rms_threshold = 0.2
      [[slits]]
          polyorder = 5
          maxshift = 0.5
          sigdetect = 50
      [[tilts]]
          tracethresh = 10, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 10
  [scienceframe]
      exprng = 20, None
      [[process]]
          satpix = nothing
          sigclip = 20.0

MAGELLAN MAGE
-------------
Alterations to the default parameters are::

  [rdx]
      spectrograph = magellan_mage
  [calibrations]
      [[biasframe]]
          useframe = overscan
      [[darkframe]]
          exprng = 20, None
      [[arcframe]]
          number = 1
          exprng = 20, None
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 3
      [[traceframe]]
          number = 3
      [[standardframe]]
          number = 1
          exprng = None, 20
      [[wavelengths]]
          echelle = True
          ech_sigrej = 3.0
          lamps = ThAr
          nonlinear_counts = 64879.65
          rms_threshold = 0.2
      [[slits]]
          polyorder = 5
          maxshift = 3.0
          pcatype = order
      [[tilts]]
          tracethresh = 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
  [scienceframe]
      exprng = 20, None
      [[process]]
          satpix = nothing
          sigclip = 20.0

KECK HIRES_R
------------
Alterations to the default parameters are::

  [rdx]
      spectrograph = keck_hires_red
  [calibrations]
      [[biasframe]]
          number = 5
      [[arcframe]]
          number = 1
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
      [[traceframe]]
          number = 3
      [[standardframe]]
          number = 1
          exprng = None, 600
      [[wavelengths]]
          echelle = True
          ech_sigrej = 3.0
          lamps = ThAr
          nonlinear_counts = 56360.1
          rms_threshold = 0.25
      [[slits]]
          polyorder = 5
          maxshift = 0.5
          sigdetect = 600.0
  [scienceframe]
      exprng = 600, None
      [[process]]
          satpix = nothing
          sigclip = 20.0

LBT MODS1R
----------
Alterations to the default parameters are::

  [rdx]
      spectrograph = lbt_mods1r
  [calibrations]
      [[biasframe]]
          number = 5
          exprng = None, 1
      [[darkframe]]
          exprng = 999999, None
      [[arcframe]]
          number = 1
          exprng = None, 60
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          exprng = 0, None
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 5
          exprng = 0, None
      [[standardframe]]
          number = 1
          exprng = 1, 200
      [[wavelengths]]
          lamps = OH_MODS
          nonlinear_counts = 64879.65
          fwhm = 10.0
          rms_threshold = 1.0
          n_first = 1
      [[slits]]
          sigdetect = 300
      [[tilts]]
          maxdev_tracefit = 0.02
          spat_order = 5
          spec_order = 5
          maxdev2d = 0.02
  [scienceframe]
      exprng = 200, None

LBT MODS1B
----------
Alterations to the default parameters are::

  [rdx]
      spectrograph = lbt_mods1b
  [calibrations]
      [[biasframe]]
          number = 5
          exprng = None, 1
      [[darkframe]]
          exprng = 999999, None
      [[arcframe]]
          number = 1
          exprng = None, 60
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          exprng = 0, None
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 5
          exprng = 0, None
      [[standardframe]]
          number = 1
          exprng = 1, 200
      [[wavelengths]]
          lamps = XeI, ArII, ArI, NeI, KrI
          nonlinear_counts = 64879.65
          rms_threshold = 0.2
          n_first = 1
      [[slits]]
          sigdetect = 300
      [[tilts]]
          maxdev_tracefit = 0.02
          spec_order = 5
          maxdev2d = 0.02
  [scienceframe]
      exprng = 200, None

LBT MODS2R
----------
Alterations to the default parameters are::

  [rdx]
      spectrograph = lbt_mods2r
  [calibrations]
      [[biasframe]]
          number = 5
          exprng = None, 1
      [[darkframe]]
          exprng = 999999, None
      [[arcframe]]
          number = 1
          exprng = None, 60
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          exprng = 0, None
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 5
          exprng = 0, None
      [[standardframe]]
          number = 1
          exprng = 1, 200
      [[wavelengths]]
          lamps = OH_MODS
          nonlinear_counts = 64879.65
          fwhm = 10.0
          rms_threshold = 1.0
          n_first = 1
      [[slits]]
          sigdetect = 300
      [[tilts]]
          maxdev_tracefit = 0.02
          spec_order = 5
          maxdev2d = 0.02
  [scienceframe]
      exprng = 200, None

LBT MODS2B
----------
Alterations to the default parameters are::

  [rdx]
      spectrograph = lbt_mods2b
  [calibrations]
      [[biasframe]]
          number = 5
          exprng = None, 1
      [[darkframe]]
          exprng = 999999, None
      [[arcframe]]
          number = 1
          exprng = None, 60
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
          exprng = 0, None
      [[pinholeframe]]
          exprng = 999999, None
      [[traceframe]]
          number = 5
          exprng = 0, None
      [[standardframe]]
          number = 1
          exprng = 1, 200
      [[wavelengths]]
          lamps = XeI, ArII, ArI, NeI, KrI
          nonlinear_counts = 64879.65
          rms_threshold = 0.2
          n_first = 1
      [[slits]]
          sigdetect = 300
      [[tilts]]
          maxdev_tracefit = 0.02
          spec_order = 5
          maxdev2d = 0.02
  [scienceframe]
      exprng = 200, None

VLT vlt_fors2
-------------
Alterations to the default parameters are::

  [rdx]
      spectrograph = vlt_fors2
  [calibrations]
      [[biasframe]]
          number = 5
      [[arcframe]]
          number = 1
          [[[process]]]
              sigrej = -1
      [[pixelflatframe]]
          number = 5
      [[traceframe]]
          number = 3
      [[standardframe]]
          number = 1
      [[flatfield]]
          illumflatten = False
          tweak_slits_thresh = 0.9
      [[wavelengths]]
          lamps = HeI, ArI
          sigdetect = 10.0
          rms_threshold = 0.25
      [[slits]]
          maxshift = 0.5
          sigdetect = 50.0
      [[tilts]]
          tracethresh = 25.0
  [flexure]
      method = boxcar

