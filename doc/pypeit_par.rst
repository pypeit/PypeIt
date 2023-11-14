.. highlight:: rest

.. include:: include/links.rst

.. _parameters:

=====================
User-level Parameters
=====================

PypeIt allows you to customize its execution without having to change the
code directly.

Although not ubiquitous, most optional arguments of PypeIt's algorithms
are contained within the :class:`~pypeit.par.pypeitpar.PypeItPar`
superset.  PypeIt uses the `configobj`_ class to parse the user-supplied
arguments  in the :ref:`pypeit_file` into an instance of
:class:`~pypeit.par.pypeitpar.PypeItPar` that is passed to all of
PypeIt's main modules.  The syntax used to set parameters using the
:ref:`pypeit_file` is important and the nesting of the parameter changes
must match the `Current PypeItPar Parameter Hierarchy`_.

.. _parameter-precedence:

Parameter Precedence
====================

The parameter changes also follow a specific precedence.  From lowest to highest,
the precedence order is as follows:

    - **Defaults**: The parameter tables below, starting with :ref:`pypeitpar`,
      all provide the *global defaults* for each parameter.

    - **Instrument-specific parameters**: Each
      :class:`~pypeit.spectrographs.spectrograph.Spectrograph` derived class
      (e.g., :class:`~pypeit.spectrographs.shane_kast.ShaneKastSpectrograph`)
      provides its own default values for
      :class:`~pypeit.par.pypeitpar.PypeItPar`, as defined by its
      ``default_pypeit_par`` method.  This allows the developers to define
      parameters as a *general expectation* for what works for each
      spectrograph.  For example, see
      :func:`~pypeit.spectrographs.shane_kast.ShaneKastSpectrograph.default_pypeit_par`
      for Shane/Kast.  All of the default changes made for each spectrograph are
      listed :ref:`here<instr_par>`.  Importantly, the parameters tabulated there
      are *not* required to be put in your :ref:`pypeit_file` if you're trying
      to reduce data from that instrument.  Only those parameters that the user
      wishes to be different from the default *as set by their specified
      instrument* need to be changed via the :ref:`pypeit_file`.  

    - **Configuration-specific parameters**: Each
      :class:`~pypeit.spectrographs.spectrograph.Spectrograph` derived class
      (e.g., :class:`~pypeit.spectrographs.shane_kast.ShaneKastSpectrograph`)
      also defines default parameters to use for specific instrument
      configurations via its ``config_specific_par`` method.  This allows the
      code to automatically define, e.g., the archived arc spectrum used for
      wavelength calibration given the grating used.  For example, see
      :func:`~pypeit.spectrographs.shane_kast.ShaneKastBlueSpectrograph.config_specific_par`
      for Shane/Kast.  These configuration-specific parameters are currently not
      documented here; however, they can be viewed by looking at the source code
      display in the API documentation.

    - **User-specified parameters**: Finally, parameters defined by the user in
      the :ref:`pypeit_file` take ultimate precedence.

.. warning::

    Default values of parameters that actually point to data files provided by
    PypeIt (e.g. the ``spectrum`` parameter for
    :class:`~pypeit.par.pypeitpar.FlexurePar`) in its root directory will point
    to the relevant location on disk of whoever generated the documentation,
    which will be different for your installation.

.. _change_par:

How to change a parameter
=========================

To change a parameter, set its value in the :ref:`parameter_block` of the
:ref:`pypeit_file`.  The *syntax* of the configuration block is important
(particularly the number of square brackets used in the parameter hierarchy),
but the indentation is not.  The indentation will just make the block easier to
read.  The :ref:`pypeit_file` :ref:`parameter_block` always includes the lines
that sets the spectrograph:

.. code-block:: ini

    [rdx]
        spectrograph = keck_deimos

The nesting of the PypeIt parameters is as illustrated in the `Current PypeItPar
Parameter Hierarchy`_ section below.  Here are a few examples of how to change
various parameters; for additional examples see the :ref:`instr_par` section.
Errors should be raised if you try to define a parameter that doesn't exist.

 * To change the threshold used for detecting slit/order edges, add:

   .. code-block:: ini

        [calibrations]
            [[slitedges]]
                edge_thresh = 100

 * To change the exposure time range used to identify an arc and
   flat-field frames and to increase the LA Cosmic sigma-clipping
   threshold for arc frames, add:

   .. code-block:: ini

        [calibrations]
            [[arcframe]]
                exprng = None,10
                [[process]]
                    sigclip = 6.
            [[pixelflatframe]]
                exprng = 11,30

.. _baseprocess:

How to change the image processing parameters for all frame types
=================================================================

To change the base-level image processing parameters that will be
applied to *all* frame types, you can use the ``baseprocess`` parameter
group.  This allows you to set these parameters once instead of having
to include lines in your PypeIt file for each frame type.  Any
frame-type-specific alterations can still be made and will overwrite the
base-level processing parameters.  For example, to change the
sigma-clipping level used by the LA Cosmic routine to default to 3.0 but
to use a value of 6.0 for arc frames, you can add the following to your
PypeIt file:

.. code-block:: ini

    [baseprocess]
        sigclip = 3.0
    [calibrations]
        [[arcframe]]
            [[[process]]]
                sigclip = 6.0

.. warning::

    Specifically for developers, note that ``baseprocess`` is only a "pseudo"
    parameter group and is not actually associated with any underlying PypeIt
    parameter class.  It is instead a flag for the code that parses the
    :ref:`pypeit_file` to distribute the associated image processing parameters
    to all of the frame-specific parameter sets.



Current PypeItPar Parameter Hierarchy
=====================================

| :ref:`pypeitpar`
|     ``[rdx]``: :ref:`reduxpar`
|     ``[calibrations]``: :ref:`calibrationspar`
|         ``[[biasframe]]``: :ref:`framegrouppar`
|             ``[[[process]]]``: :ref:`processimagespar`
|         ``[[darkframe]]``: :ref:`framegrouppar`
|             ``[[[process]]]``: :ref:`processimagespar`
|         ``[[arcframe]]``: :ref:`framegrouppar`
|             ``[[[process]]]``: :ref:`processimagespar`
|         ``[[tiltframe]]``: :ref:`framegrouppar`
|             ``[[[process]]]``: :ref:`processimagespar`
|         ``[[pixelflatframe]]``: :ref:`framegrouppar`
|             ``[[[process]]]``: :ref:`processimagespar`
|         ``[[pinholeframe]]``: :ref:`framegrouppar`
|             ``[[[process]]]``: :ref:`processimagespar`
|         ``[[alignframe]]``: :ref:`framegrouppar`
|             ``[[[process]]]``: :ref:`processimagespar`
|         ``[[alignment]]``: :ref:`alignpar`
|         ``[[traceframe]]``: :ref:`framegrouppar`
|             ``[[[process]]]``: :ref:`processimagespar`
|         ``[[illumflatframe]]``: :ref:`framegrouppar`
|             ``[[[process]]]``: :ref:`processimagespar`
|         ``[[lampoffflatsframe]]``: :ref:`framegrouppar`
|             ``[[[process]]]``: :ref:`processimagespar`
|         ``[[skyframe]]``: :ref:`framegrouppar`
|             ``[[[process]]]``: :ref:`processimagespar`
|         ``[[standardframe]]``: :ref:`framegrouppar`
|             ``[[[process]]]``: :ref:`processimagespar`
|         ``[[flatfield]]``: :ref:`flatfieldpar`
|         ``[[wavelengths]]``: :ref:`wavelengthsolutionpar`
|         ``[[slitedges]]``: :ref:`edgetracepar`
|         ``[[tilts]]``: :ref:`wavetiltspar`
|     ``[scienceframe]``: :ref:`framegrouppar`
|         ``[[process]]``: :ref:`processimagespar`
|     ``[reduce]``: :ref:`reducepar`
|         ``[[findobj]]``: :ref:`findobjpar`
|         ``[[skysub]]``: :ref:`skysubpar`
|         ``[[extraction]]``: :ref:`extractionpar`
|         ``[[cube]]``: :ref:`cubepar`
|         ``[[slitmask]]``: :ref:`slitmaskpar`
|     ``[flexure]``: :ref:`flexurepar`
|     ``[fluxcalib]``: :ref:`fluxcalibratepar`
|     ``[coadd1d]``: :ref:`coadd1dpar`
|     ``[coadd2d]``: :ref:`coadd2dpar`
|     ``[sensfunc]``: :ref:`sensfuncpar`
|         ``[[UVIS]]``: :ref:`sensfuncuvispar`
|         ``[[IR]]``: :ref:`telluricpar`
|     ``[telluric]``: :ref:`telluricpar`
|     ``[collate1d]``: :ref:`collate1dpar`

----

.. _pypeitpar:

PypeItPar Keywords
------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.PypeItPar`

================  ===============================================  =======  ============================  ======================================================================================================================================================================================================================================================================================
Key               Type                                             Options  Default                       Description                                                                                                                                                                                                                                                                           
================  ===============================================  =======  ============================  ======================================================================================================================================================================================================================================================================================
``calibrations``  :class:`~pypeit.par.pypeitpar.CalibrationsPar`   ..       `CalibrationsPar Keywords`_   Parameters for the calibration algorithms                                                                                                                                                                                                                                             
``coadd1d``       :class:`~pypeit.par.pypeitpar.Coadd1DPar`        ..       `Coadd1DPar Keywords`_        Par set to control 1D coadds.  Only used in the after-burner script.                                                                                                                                                                                                                  
``coadd2d``       :class:`~pypeit.par.pypeitpar.Coadd2DPar`        ..       `Coadd2DPar Keywords`_        Par set to control 2D coadds.  Only used in the after-burner script.                                                                                                                                                                                                                  
``collate1d``     :class:`~pypeit.par.pypeitpar.Collate1DPar`      ..       `Collate1DPar Keywords`_      Par set to control collating 1d spectra.  Only used in the after-burner script.                                                                                                                                                                                                       
``flexure``       :class:`~pypeit.par.pypeitpar.FlexurePar`        ..       `FlexurePar Keywords`_        Parameters used by the flexure-correction procedure.  Flexure corrections are not performed by default.  To turn on, either set the parameters in the 'flexure' parameter group or set 'flexure = True' in the 'rdx' parameter group to use the default flexure-correction parameters.
``fluxcalib``     :class:`~pypeit.par.pypeitpar.FluxCalibratePar`  ..       `FluxCalibratePar Keywords`_  Parameters used by the flux-calibration procedure.  Flux calibration is not performed by default.  To turn on, either set the parameters in the 'fluxcalib' parameter group or set 'fluxcalib = True' in the 'rdx' parameter group to use the default flux-calibration parameters.    
``rdx``           :class:`~pypeit.par.pypeitpar.ReduxPar`          ..       `ReduxPar Keywords`_          PypeIt reduction rules.                                                                                                                                                                                                                                                               
``reduce``        :class:`~pypeit.par.pypeitpar.ReducePar`         ..       `ReducePar Keywords`_         Parameters determining sky-subtraction, object finding, and extraction                                                                                                                                                                                                                
``scienceframe``  :class:`~pypeit.par.pypeitpar.FrameGroupPar`     ..       `FrameGroupPar Keywords`_     The frames and combination rules for the science observations                                                                                                                                                                                                                         
``sensfunc``      :class:`~pypeit.par.pypeitpar.SensFuncPar`       ..       `SensFuncPar Keywords`_       Par set to control sensitivity function computation.  Only used in the after-burner script.                                                                                                                                                                                           
``telluric``      :class:`~pypeit.par.pypeitpar.TelluricPar`       ..       `TelluricPar Keywords`_       Par set to control telluric fitting.  Only used in the pypeit_sensfunc and pypeit_telluric after-burner scripts.                                                                                                                                                                      
================  ===============================================  =======  ============================  ======================================================================================================================================================================================================================================================================================


----

.. _calibrationspar:

CalibrationsPar Keywords
------------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.CalibrationsPar`

=====================  ====================================================  =======  =================================  =================================================================================================================================================================================================================================================
Key                    Type                                                  Options  Default                            Description                                                                                                                                                                                                                                      
=====================  ====================================================  =======  =================================  =================================================================================================================================================================================================================================================
``alignframe``         :class:`~pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the align frames                                                                                                                                                                                            
``alignment``          :class:`~pypeit.par.pypeitpar.AlignPar`               ..       `AlignPar Keywords`_               Define the procedure for the alignment of traces                                                                                                                                                                                                 
``arcframe``           :class:`~pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the wavelength calibration                                                                                                                                                                                  
``biasframe``          :class:`~pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the bias correction                                                                                                                                                                                         
``bpm_usebias``        bool                                                  ..       False                              Make a bad pixel mask from bias frames? Bias frames must be provided.                                                                                                                                                                            
``calib_dir``          str                                                   ..       ``Calibrations``                   The name of the directory for the processed calibration frames.  The host path for the directory is set by the redux_path (see :class:`~pypeit.par.pypeitpar.ReduxPar`).  Beware that success when changing the default value is not well tested!
``darkframe``          :class:`~pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the dark-current correction                                                                                                                                                                                 
``flatfield``          :class:`~pypeit.par.pypeitpar.FlatFieldPar`           ..       `FlatFieldPar Keywords`_           Parameters used to set the flat-field procedure                                                                                                                                                                                                  
``illumflatframe``     :class:`~pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the illumination flat                                                                                                                                                                                       
``lampoffflatsframe``  :class:`~pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the lamp off flats                                                                                                                                                                                          
``pinholeframe``       :class:`~pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the pinholes                                                                                                                                                                                                
``pixelflatframe``     :class:`~pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the pixel flat                                                                                                                                                                                              
``raise_chk_error``    bool                                                  ..       True                               Raise an error if the calibration check fails                                                                                                                                                                                                    
``skyframe``           :class:`~pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the sky background observations                                                                                                                                                                             
``slitedges``          :class:`~pypeit.par.pypeitpar.EdgeTracePar`           ..       `EdgeTracePar Keywords`_           Slit-edge tracing parameters                                                                                                                                                                                                                     
``standardframe``      :class:`~pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the spectrophotometric standard observations                                                                                                                                                                
``tiltframe``          :class:`~pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for the wavelength tilts                                                                                                                                                                                        
``tilts``              :class:`~pypeit.par.pypeitpar.WaveTiltsPar`           ..       `WaveTiltsPar Keywords`_           Define how to trace the slit tilts using the trace frames                                                                                                                                                                                        
``traceframe``         :class:`~pypeit.par.pypeitpar.FrameGroupPar`          ..       `FrameGroupPar Keywords`_          The frames and combination rules for images used for slit tracing                                                                                                                                                                                
``wavelengths``        :class:`~pypeit.par.pypeitpar.WavelengthSolutionPar`  ..       `WavelengthSolutionPar Keywords`_  Parameters used to derive the wavelength solution                                                                                                                                                                                                
=====================  ====================================================  =======  =================================  =================================================================================================================================================================================================================================================


----

.. _alignpar:

AlignPar Keywords
-----------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.AlignPar`

===============  =============  =======  ========  ================================================================================================================================================================================================================================================================================
Key              Type           Options  Default   Description                                                                                                                                                                                                                                                                     
===============  =============  =======  ========  ================================================================================================================================================================================================================================================================================
``locations``    list, ndarray  ..       0.0, 1.0  Locations of the bars, in a list, specified as a fraction of the slit width                                                                                                                                                                                                     
``snr_thresh``   int, float     ..       1.0       S/N ratio threshold for finding an alignment trace. This should be a low number to ensure that the algorithm finds all bars. The algorithm will then only use the N most significant detections, where N is the number of elements specified in the "locations" keyword argument
``trace_npoly``  int            ..       4         Order of the polynomial to use when fitting the trace of a single bar                                                                                                                                                                                                           
``trim_edge``    list           ..       0, 0      Trim the slit by this number of pixels left/right before finding alignment bars                                                                                                                                                                                                 
===============  =============  =======  ========  ================================================================================================================================================================================================================================================================================


----

.. _flatfieldpar:

FlatFieldPar Keywords
---------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.FlatFieldPar`

==========================  =================  =================================  ===========  ================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                         Type               Options                            Default      Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
==========================  =================  =================================  ===========  ================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``fit_2d_det_response``     bool               ..                                 False        Set this variable to True if you want to compute and account for the detector response in the flatfield image. Note that ``detector response`` refers to pixel sensitivity variations that primarily depend on (x,y) detector coordinates. In most cases, the default 2D bspline is sufficient to account for detector response (i.e. set this parameter to False). Note that this correction will _only_ be performed for the spectrographs that have a dedicated response correction implemented. Currently,this correction is only implemented for Keck+KCWI.
``illum_iter``              int                ..                                 0            The number of rejection iterations to perform when constructing the slit-illumination profile.  No rejection iterations are performed if 0.  WARNING: Functionality still being tested.                                                                                                                                                                                                                                                                                                                                                                         
``illum_rej``               int, float         ..                                 5.0          The sigma threshold used in the rejection iterations used to refine the slit-illumination profile.  Rejection iterations are only performed if ``illum_iter > 0``.                                                                                                                                                                                                                                                                                                                                                                                              
``method``                  str                ``bspline``, ``skip``              ``bspline``  Method used to flat field the data; use skip to skip flat-fielding.  Options are: None, bspline, skip                                                                                                                                                                                                                                                                                                                                                                                                                                                           
``pixelflat_file``          str                ..                                 ..           Filename of the image to use for pixel-level field flattening                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``pixelflat_max_wave``      int, float         ..                                 ..           All values of the normalized pixel flat are set to 1 for wavelengths above this value.                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
``pixelflat_min_wave``      int, float         ..                                 ..           All values of the normalized pixel flat are set to 1 for wavelengths below this value.                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
``rej_sticky``              bool               ..                                 False        Propagate the rejected pixels through the stages of the flat-field fitting (i.e, from the spectral fit, to the spatial fit, and finally to the 2D residual fit).  If False, pixels rejected in each stage are included in each subsequent stage.                                                                                                                                                                                                                                                                                                                
``saturated_slits``         str                ``crash``, ``mask``, ``continue``  ``crash``    Behavior when a slit is encountered with a large fraction of saturated pixels in the flat-field.  The options are: 'crash' - Raise an error and halt the data reduction; 'mask' - Mask the slit, meaning no science data will be extracted from the slit; 'continue' - ignore the flat-field correction, but continue with the reduction.                                                                                                                                                                                                                       
``slit_illum_finecorr``     bool               ..                                 True         If True, a fine correction to the spatial illumination profile will be performed. The fine correction is a low order 2D polynomial fit to account for a gradual change to the spatial illumination profile as a function of wavelength.                                                                                                                                                                                                                                                                                                                         
``slit_illum_pad``          int, float         ..                                 5.0          The number of pixels to pad the slit edges when constructing the slit-illumination profile. Single value applied to both edges.                                                                                                                                                                                                                                                                                                                                                                                                                                 
``slit_illum_ref_idx``      int                ..                                 0            The index of a reference slit (0-indexed) used for estimating the relative spectral sensitivity (or the relative blaze). This parameter is only used if ``slit_illum_relative = True``.                                                                                                                                                                                                                                                                                                                                                                         
``slit_illum_relative``     bool               ..                                 False        Generate an image of the relative spectral illumination for a multi-slit setup.  If you set ``use_slitillum = True`` for any of the frames that use the flatfield model, this *must* be set to True. Currently, this is only used for IFU reductions.                                                                                                                                                                                                                                                                                                           
``slit_illum_smooth_npix``  int                ..                                 10           The number of pixels used to determine smoothly varying relative weights is given by ``nspec/slit_illum_smooth_npix``, where nspec is the number of spectral pixels.                                                                                                                                                                                                                                                                                                                                                                                            
``slit_trim``               int, float, tuple  ..                                 3.0          The number of pixels to trim each side of the slit when selecting pixels to use for fitting the spectral response function.  Single values are used for both slit edges; a two-tuple can be used to trim the left and right sides differently.                                                                                                                                                                                                                                                                                                                  
``spat_samp``               int, float         ..                                 5.0          Spatial sampling for slit illumination function. This is the width of the median filter in pixels used to determine the slit illumination function, and thus sets the minimum scale on which the illumination function will have features.                                                                                                                                                                                                                                                                                                                      
``spec_samp_coarse``        int, float         ..                                 50.0         bspline break point spacing in units of pixels for 2-d bspline-polynomial fit to flat field image residuals. This should be a large number unless you are trying to fit a sky flat with lots of narrow spectral features.                                                                                                                                                                                                                                                                                                                                       
``spec_samp_fine``          int, float         ..                                 1.2          bspline break point spacing in units of pixels for spectral fit to flat field blaze function.                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``tweak_slits``             bool               ..                                 True         Use the illumination flat field to tweak the slit edges. This will work even if illumflatten is set to False                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``tweak_slits_maxfrac``     float              ..                                 0.1          If tweak_slit is True, this sets the maximum fractional amount (of a slits width) allowed for trimming each (i.e. left and right) slit boundary, i.e. the default is 10% which means slits would shrink or grow by at most 20% (10% on each side)                                                                                                                                                                                                                                                                                                               
``tweak_slits_thresh``      float              ..                                 0.93         If tweak_slits is True, this sets the illumination function threshold used to tweak the slit boundaries based on the illumination flat. It should be a number less than 1.0                                                                                                                                                                                                                                                                                                                                                                                     
``twod_fit_npoly``          int                ..                                 ..           Order of polynomial used in the 2D bspline-polynomial fit to flat-field image residuals. The code determines the order of these polynomials to each slit automatically depending on the slit width, which is why the default is None. Alter this paramter at your own risk!                                                                                                                                                                                                                                                                                     
==========================  =================  =================================  ===========  ================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _edgetracepar:

EdgeTracePar Keywords
---------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.EdgeTracePar`

===========================  ================  ===========================================  ==============  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                          Type              Options                                      Default         Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
===========================  ================  ===========================================  ==============  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``add_missed_orders``        bool              ..                                           False           If orders are not detected by the automated edge tracing, attempt to add them based on their expected positions on on the detector.  Echelle spectrographs only.                                                                                                                                                                                                                                                                                                                                                                                                                                      
``add_predict``              str               ..                                           ``nearest``     Sets the method used to predict the shape of the left and right traces for a user-defined slit inserted.  Options are (1) ``straight`` inserts traces with a constant spatial pixels position, (2) ``nearest`` inserts traces with a form identical to the automatically identified trace at the nearest spatial position to the inserted slit, or (3) ``pca`` uses the PCA decomposition to predict the shape of the traces.                                                                                                                                                                         
``add_slits``                str, list         ..                                           ..              Add one or more user-defined slits.  The syntax to define a slit to add is: 'det:spec:spat_left:spat_right' where det=detector, spec=spectral pixel, spat_left=spatial pixel of left slit boundary, and spat_righ=spatial pixel of right slit boundary.  For example, '2:2000:2121:2322,3:2000:1201:1500' will add a slit to detector 2 passing through spec=2000 extending spatially from 2121 to 2322 and another on detector 3 at spec=2000 extending from 1201 to 1500.                                                                                                                           
``auto_pca``                 bool              ..                                           True            During automated tracing, attempt to construct a PCA decomposition of the traces. When True, the edge traces resulting from the initial detection, centroid refinement, and polynomial fitting must meet a set of criteria for performing the pca; see :func:`pypeit.edgetrace.EdgeTraceSet.can_pca`.  If False, the ``sync_predict`` parameter *cannot* be set to ``pca``; if it is not, the value is set to ``nearest`` and a warning is issued when validating the parameter set.                                                                                                                  
``bound_detector``           bool              ..                                           False           When the code is ready to synchronize the left/right trace edges, the traces should have been constructed, vetted, and cleaned. This can sometimes lead to *no* valid traces. This parameter dictates what to do next. If ``bound_detector`` is True, the code will artificially add left and right edges that bound the detector; if False, the code identifies the slit-edge tracing as being unsuccessful, warns the user, and ends gracefully. Note that setting ``bound_detector`` to True is needed for some long-slit data where the slit edges are, in fact, beyond the edges of the detector.
``clip``                     bool              ..                                           True            Remove traces flagged as bad, instead of only masking them.  This is currently only used by :func:`~pypeit.edgetrace.EdgeTraceSet.centroid_refine`.                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``det_buffer``               int               ..                                           5               The minimum separation between the detector edges and a slit edge for any added edge traces.  Must be positive.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
``det_min_spec_length``      int, float        ..                                           0.33            The minimum spectral length (as a fraction of the detector size) of a trace determined by direct measurements of the detector data (as opposed to what should be included in any modeling approach; see fit_min_spec_length).                                                                                                                                                                                                                                                                                                                                                                         
``dlength_range``            int, float        ..                                           ..              Similar to ``minimum_slit_dlength``, but constrains the *fractional* change in the slit length as a function of wavelength.  For example, a value of 0.2 means that slit length should not vary more than 20%as a function of wavelength.                                                                                                                                                                                                                                                                                                                                                             
``edge_detect_clip``         int, float        ..                                           ..              Sigma clipping level for peaks detected in the collapsed, Sobel-filtered significance image.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
``edge_thresh``              int, float        ..                                           20.0            Threshold for finding edges in the Sobel-filtered significance image.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
``exclude_regions``          list, str         ..                                           ..              User-defined regions to exclude from the slit tracing. To set this parameter, the text should be a comma separated list of pixel ranges (in the x direction) to be excluded and the detector number. For example, the following string 1:0:20,1:300:400  would select two regions in det=1 between pixels 0 and 20 and between 300 and 400.                                                                                                                                                                                                                                                           
``filt_iter``                int               ..                                           0               Number of median-filtering iterations to perform on sqrt(trace) image before applying to Sobel filter to detect slit/order edges.                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``fit_function``             str               ``polynomial``, ``legendre``, ``chebyshev``  ``legendre``    Function fit to edge measurements.  Options are: polynomial, legendre, chebyshev                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
``fit_maxdev``               int, float        ..                                           5.0             Maximum deviation between the fitted and measured edge position for rejection in spatial pixels.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
``fit_maxiter``              int               ..                                           25              Maximum number of rejection iterations during edge fitting.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
``fit_min_spec_length``      float             ..                                           0.6             Minimum unmasked spectral length of a traced slit edge to use in any modeling procedure (polynomial fitting or PCA decomposition).                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``fit_niter``                int               ..                                           1               Number of iterations of re-measuring and re-fitting the edge data; see :func:`~pypeit.core.trace.fit_trace`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
``fit_order``                int               ..                                           5               Order of the function fit to edge measurements.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
``follow_span``              int               ..                                           20              In the initial connection of spectrally adjacent edge detections, this sets the number of previous spectral rows to consider when following slits forward.                                                                                                                                                                                                                                                                                                                                                                                                                                            
``fwhm_gaussian``            int, float        ..                                           3.0             The `fwhm` parameter to use when using Gaussian weighting in :func:`~pypeit.core.trace.fit_trace` when refining the PCA predictions of edges.  See description :func:`~pypeit.core.trace.peak_trace`.                                                                                                                                                                                                                                                                                                                                                                                                 
``fwhm_uniform``             int, float        ..                                           3.0             The `fwhm` parameter to use when using uniform weighting in :func:`~pypeit.core.trace.fit_trace` when refining the PCA predictions of edges.  See description of :func:`~pypeit.core.trace.peak_trace`.                                                                                                                                                                                                                                                                                                                                                                                               
``gap_offset``               int, float        ..                                           5.0             Offset (pixels) used for the slit edge gap width when inserting slit edges (see `sync_center`) or when nudging predicted slit edges to avoid slit overlaps.  This should be larger than `minimum_slit_gap` when converted to arcseconds.                                                                                                                                                                                                                                                                                                                                                              
``left_right_pca``           bool              ..                                           False           Construct a PCA decomposition for the left and right traces separately.  This can be important for cross-dispersed echelle spectrographs (e.g., Keck-NIRES)                                                                                                                                                                                                                                                                                                                                                                                                                                           
``length_range``             int, float        ..                                           ..              Allowed range in slit length compared to the median slit length.  For example, a value of 0.3 means that slit lengths should not vary more than 30%.  Relatively shorter or longer slits are masked or clipped.  Most useful for echelle or multi-slit data where the slits should have similar or identical lengths.                                                                                                                                                                                                                                                                                 
``maskdesign_filename``      str, list         ..                                           ..              Mask design info contained in this file or files (comma separated)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``maskdesign_maxsep``        int, float        ..                                           50              Maximum allowed offset in pixels between the slit edges defined by the slit-mask design and the traced edges.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
``maskdesign_sigrej``        int, float        ..                                           3               Number of sigma for sigma-clipping rejection during slit-mask design matching.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``maskdesign_step``          int, float        ..                                           1               Step in pixels used to generate a list of possible offsets (within +/- `maskdesign_maxsep`) between the slit edges defined by the mask design and the traced edges.                                                                                                                                                                                                                                                                                                                                                                                                                                   
``match_tol``                int, float        ..                                           3.0             Same-side slit edges below this separation in pixels are considered part of the same edge.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
``max_nudge``                int, float        ..                                           ..              If parts of any (predicted) trace fall off the detector edge, allow them to be nudged away from the detector edge up to and including this maximum number of pixels.  If None, no limit is set; otherwise should be 0 or larger.                                                                                                                                                                                                                                                                                                                                                                      
``max_shift_abs``            int, float        ..                                           0.5             Maximum spatial shift in pixels between an input edge location and the recentroided value.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
``max_shift_adj``            int, float        ..                                           0.15            Maximum spatial shift in pixels between the edges in adjacent spectral positions.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``max_spat_error``           int, float        ..                                           ..              Maximum error in the spatial position of edges in pixels.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``minimum_slit_dlength``     int, float        ..                                           ..              Minimum *change* in the slit length (arcsec) as a function of wavelength in arcsec.  This is mostly meant to catch cases when the polynomial fit to the detected edges becomes ill-conditioned (e.g., when the slits run off the edge of the detector) and leads to wild traces.  If reducing the order of the polynomial (``fit_order``) does not help, try using this to remove poorly constrained slits.                                                                                                                                                                                           
``minimum_slit_gap``         int, float        ..                                           ..              Minimum slit gap in arcsec.  Gaps between slits are determined by the median difference between the right and left edge locations of adjacent slits.  Slits with small gaps are merged by removing the intervening traces.If None, no minimum slit gap is applied.  This should be smaller than `gap_offset` when converted to pixels.                                                                                                                                                                                                                                                                
``minimum_slit_length``      int, float        ..                                           ..              Minimum slit length in arcsec.  Slit lengths are determined by the median difference between the left and right edge locations for the unmasked trace locations.  This is used to identify traces that are *erroneously* matched together to form slits.  Short slits are expected to be ignored or removed (see  ``clip``).  If None, no minimum slit length applied.                                                                                                                                                                                                                                
``minimum_slit_length_sci``  int, float        ..                                           ..              Minimum slit length in arcsec for a science slit.  Slit lengths are determined by the median difference between the left and right edge locations for the unmasked trace locations.  Used in combination with ``minimum_slit_length``, this parameter is used to identify box or alignment slits; i.e., those slits that are shorter than ``minimum_slit_length_sci`` but larger than ``minimum_slit_length`` are box/alignment slits.  Box slits are *never* removed (see ``clip``), but no spectra are extracted from them.  If None, no minimum science slit length is applied.                    
``niter_gaussian``           int               ..                                           6               The number of iterations of :func:`~pypeit.core.trace.fit_trace` to use when using Gaussian weighting.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
``niter_uniform``            int               ..                                           9               The number of iterations of :func:`~pypeit.core.trace.fit_trace` to use when using uniform weighting.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
``order_match``              int, float        ..                                           ..              For echelle spectrographs, this is the tolerance allowed for matching identified "slits" to echelle orders. Must be in the fraction of the detector spatial scale (i.e., a value of 0.05 means that the order locations must be within 5% of the expected value).  If None, no limit is used.                                                                                                                                                                                                                                                                                                         
``order_offset``             int, float        ..                                           ..              Offset to introduce to the expected order positions to improve the match for this specific data. This is an additive offset to the measured slit positions; i.e., this should minimize the difference between the expected order positions and ``self.slit_spatial_center() + offset``. Must be in the fraction of the detector spatial scale. If None, no offset is applied.                                                                                                                                                                                                                         
``overlap``                  bool              ..                                           False           Assume slits identified as abnormally short are actually due to overlaps between adjacent slits/orders.  If set to True, you *must* have also used ``length_range`` to identify left-right edge pairs that have an abnormally short separation.  For those short slits, the code attempts to convert the short slits into slit gaps.  This is particularly useful for blue orders in Keck-HIRES data.                                                                                                                                                                                                 
``pad``                      int               ..                                           0               Integer number of pixels to consider beyond the slit edges when selecting pixels that are 'on' the slit.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
``pca_function``             str               ``polynomial``, ``legendre``, ``chebyshev``  ``polynomial``  Type of function fit to the PCA coefficients for each component.  Options are: polynomial, legendre, chebyshev                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``pca_maxiter``              int               ..                                           25              Maximum number of rejection iterations when fitting the PCA coefficients.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``pca_maxrej``               int               ..                                           1               Maximum number of PCA coefficients rejected during a given fit iteration.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``pca_min_edges``            int               ..                                           4               Minimum number of edge traces required to perform a PCA decomposition of the trace form.  If left_right_pca is True, this minimum applies to the number of left and right traces separately.                                                                                                                                                                                                                                                                                                                                                                                                          
``pca_n``                    int               ..                                           ..              The number of PCA components to keep, which must be less than the number of detected traces.  If not provided, determined by calculating the minimum number of components required to explain a given percentage of variance in the edge data; see `pca_var_percent`.                                                                                                                                                                                                                                                                                                                                 
``pca_order``                int               ..                                           2               Order of the function fit to the PCA coefficients.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``pca_sigrej``               int, float, list  ..                                           2.0, 2.0        Sigma rejection threshold for fitting PCA components. Individual numbers are used for both lower and upper rejection. A list of two numbers sets these explicitly (e.g., [2., 3.]).                                                                                                                                                                                                                                                                                                                                                                                                                   
``pca_var_percent``          int, float        ..                                           99.8            The percentage (i.e., not the fraction) of the variance in the edge data accounted for by the PCA used to truncate the number of PCA coefficients to keep (see `pca_n`).  Ignored if `pca_n` is provided directly.                                                                                                                                                                                                                                                                                                                                                                                    
``rm_slits``                 str, list         ..                                           ..              Remove one or more user-specified slits.  The syntax used to define a slit to remove is: 'det:spec:spat' where det=detector, spec=spectral pixel, spat=spatial pixel.  For example, '2:2000:2121,3:2000:1500' will remove the slit on detector 2 that contains pixel (spat,spec)=(2000,2121) and on detector 3 that contains pixel (2000,2121).                                                                                                                                                                                                                                                       
``smash_range``              list              ..                                           0.0, 1.0        Range of the slit in the spectral direction (in fractional units) to smash when searching for slit edges.  If the spectrum covers only a portion of the image, use that range.                                                                                                                                                                                                                                                                                                                                                                                                                        
``sobel_enhance``            int               ..                                           0               Enhance the sobel filtering? A value of 0 will not enhance the sobel filtering. Any other value > 0 will sum the sobel values. For example, a value of 3 will combine the sobel values for the 3 nearest pixels. This is useful when a slit edge is poorly defined (e.g. vignetted).                                                                                                                                                                                                                                                                                                                  
``sobel_mode``               str               ``nearest``, ``constant``                    ``nearest``     Mode for Sobel filtering.  Default is 'nearest'; note we find'constant' works best for DEIMOS.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``sync_center``              str               ``median``, ``nearest``, ``gap``             ``median``      Mode to use for determining the location of traces to insert.  Use `median` to use the median of the matched left and right edge pairs, `nearest` to use the length of the nearest slit, or `gap` to offset by a fixed gap width from the next slit edge.                                                                                                                                                                                                                                                                                                                                             
``sync_predict``             str               ``pca``, ``nearest``, ``auto``               ``pca``         Mode to use when predicting the form of the trace to insert.  Use `pca` to use the PCA decomposition, `nearest` to reproduce the shape of the nearest trace, or `auto` to let PypeIt decide which mode to use between `pca` and `nearest`. In general, it will first try `pca`, and if that is not possible, it will use `nearest`.                                                                                                                                                                                                                                                                   
``sync_to_edge``             bool              ..                                           True            If adding a first left edge or a last right edge, ignore `center_mode` for these edges and place them at the edge of the detector (with the relevant shape).                                                                                                                                                                                                                                                                                                                                                                                                                                          
``trace_median_frac``        int, float        ..                                           ..              After detection of peaks in the rectified Sobel-filtered image and before refitting the edge traces, the rectified image is median filtered with a kernel width of `trace_median_frac*nspec` along the spectral dimension.                                                                                                                                                                                                                                                                                                                                                                            
``trace_thresh``             int, float        ..                                           ..              After rectification and median filtering of the Sobel-filtered image (see `trace_median_frac`), values in the median-filtered image *below* this threshold are masked in the refitting of the edge trace data.  If None, no masking applied.                                                                                                                                                                                                                                                                                                                                                          
``use_maskdesign``           bool              ..                                           False           Use slit-mask designs to identify slits.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
===========================  ================  ===========================================  ==============  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _wavetiltspar:

WaveTiltsPar Keywords
---------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.WaveTiltsPar`

===================  =========================  =======  ==============  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                  Type                       Options  Default         Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
===================  =========================  =======  ==============  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``cont_rej``         int, float, list, ndarray  ..       3, 1.5          The sigma threshold for rejection.  Can be a single number or two numbers that give the low and high sigma rejection, respectively.                                                                                                                                                                                                                                                                                                                                                                                                                                          
``func2d``           str                        ..       ``legendre2d``  Type of function for 2D fit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``idsonly``          bool                       ..       False           Only use the arc lines that have an identified wavelength to trace tilts (CURRENTLY NOT USED!)                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
``maxdev2d``         int, float                 ..       0.25            Maximum absolute deviation (in units of fwhm) rejection threshold used to determines which pixels in global 2d fits to arc line tilts are rejected because they deviate from the model by more than this value                                                                                                                                                                                                                                                                                                                                                               
``maxdev_tracefit``  int, float                 ..       0.2             Maximum absolute deviation (in units of fwhm) for the legendre polynomial fits to individual arc line tilt fits during iterative trace fitting (flux weighted, then gaussian weighted)                                                                                                                                                                                                                                                                                                                                                                                       
``minmax_extrap``    list, ndarray              ..       150.0, 1000.0   Sets how far below the last measured tilt line is extrapolated in tracewave.fit_tilts()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
``nfwhm_neigh``      int, float                 ..       3.0             Required separation between neighboring arc lines for them to be considered for tilt tracing in units of the the spectral fwhm (see wavelength parset where fwhm is defined)                                                                                                                                                                                                                                                                                                                                                                                                 
``rm_continuum``     bool                       ..       False           Before tracing the line center at each spatial position, remove any low-order continuum in the 2D spectra.                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``sig_neigh``        int, float                 ..       10.0            Significance threshold for arcs to be used in line identification for the purpose of identifying neighboring lines. The tracethresh parameter above determines the significance threshold of lines that will be traced, but these lines  must be at least nfwhm_neigh fwhm away from neighboring lines. This parameter determines the significance above which  a line must be to be considered a possible colliding neighbor. A low value of sig_neigh will result in an overall  larger number of lines, which will result in more lines above tracethresh getting rejected
``sigrej2d``         int, float                 ..       3.0             Outlier rejection significance determining which pixels on a fit to an arc line tilt are rejected by the global 2D fit                                                                                                                                                                                                                                                                                                                                                                                                                                                       
``sigrej_trace``     int, float                 ..       3.0             Outlier rejection significance to determine which traced arc lines should be included in the global fit                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
``spat_order``       int, float, list, ndarray  ..       3               Order of the legendre polynomial to be fit to the tilt of an arc line. This parameter determines both the order of the *individual* arc line tilts, as well as the order of the spatial direction of the 2d legendre polynomial (spatial, spectral) that is fit to obtain a global solution for the tilts across the slit/order. This can be a single number or a list/array providing the value for each slit                                                                                                                                                               
``spec_order``       int, float, list, ndarray  ..       4               Order of the spectral direction of the 2d legendre polynomial (spatial, spectral) that is fit to obtain a global solution for the tilts across the slit/order. This can be a single number or a list/array providing the value for each slit                                                                                                                                                                                                                                                                                                                                 
``tracethresh``      int, float, list, ndarray  ..       20.0            Significance threshold for arcs to be used in tracing wavelength tilts. This can be a single number or a list/array providing the value for each slit/order.                                                                                                                                                                                                                                                                                                                                                                                                                 
===================  =========================  =======  ==============  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _wavelengthsolutionpar:

WavelengthSolutionPar Keywords
------------------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.WavelengthSolutionPar`

====================  =========================  ============================================================================  ================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                   Type                       Options                                                                       Default           Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
====================  =========================  ============================================================================  ================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``cc_local_thresh``   float                      ..                                                                            0.7               Threshold for the *local* cross-correlation coefficient, evaluated at each reidentified line,  between an input spectrum and the shifted and stretched archive spectrum above which a line must be to be considered a good line for reidentification. The local cross-correlation is evaluated at each candidate reidentified line (using a window of nlocal_cc), and is then used to score the the reidentified lines to arrive at the final set of good reidentifications.                                                                                                                                                                                                                                                                                                                                 
``cc_thresh``         float, list, ndarray       ..                                                                            0.7               Threshold for the *global* cross-correlation coefficient between an input spectrum and member of the archive required to attempt reidentification.  Spectra from the archive with a lower cross-correlation are not used for reidentification. This can be a single number or a list/array providing the value for each slit.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
``ech_norder_coeff``  int                        ..                                                                            4                 For echelle spectrographs, this is the order of the final 2d fit to the order dimension.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``ech_nspec_coeff``   int                        ..                                                                            4                 For echelle spectrographs, this is the order of the final 2d fit to the spectral dimension.  You should choose this to be the n_final of the fits to the individual orders.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``ech_separate_2d``   bool                       ..                                                                            False             For echelle spectrographs, fit the 2D solutions on separate detectors separately                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``ech_sigrej``        int, float                 ..                                                                            2.0               For echelle spectrographs, this is the sigma-clipping rejection threshold in the 2d fit to spectral and order dimensions                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``echelle``           bool                       ..                                                                            False             Is this an echelle spectrograph? If yes an additional 2-d fit wavelength fit will be performed as a function of spectral pixel and order number to improve the wavelength solution                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
``func``              str                        ..                                                                            ``legendre``      Function used for wavelength solution fits                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``fwhm``              int, float                 ..                                                                            4.0               Spectral sampling of the arc lines. This is the FWHM of an arcline in binned pixels of the input arc image                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``fwhm_fromlines``    bool                       ..                                                                            False             Estimate spectral resolution in each slit using the arc lines. If True, the estimated FWHM will override ``fwhm`` only in the determination of the wavelength solution (`i.e.`, not in WaveTilts).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
``fwhm_spat_order``   int                        ..                                                                            0                 This parameter determines the spatial polynomial order to use in the 2D polynomial fit to the FWHM of the arc lines. See also, fwhm_spec_order.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
``fwhm_spec_order``   int                        ..                                                                            1                 This parameter determines the spectral polynomial order to use in the 2D polynomial fit to the FWHM of the arc lines. See also, fwhm_spat_order.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``lamps``             list                       ..                                                                            ..                Name of one or more ions used for the wavelength calibration.  Use ``None`` for no calibration. Choose ``use_header`` to use the list of lamps recorded in the header of the arc frames (this is currently available only for Keck DEIMOS and LDT DeVeny).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``match_toler``       float                      ..                                                                            2.0               Matching tolerance in pixels when searching for new lines. This is the difference in pixels between the wavlength assigned to an arc line by an iteration of the wavelength solution to the wavelength in the line list.  This parameter is also used as the matching tolerance in pixels for a line reidentification.  A good line match must match within this tolerance to the shifted and stretched archive spectrum, and the archive wavelength solution at this match must be within match_toler dispersion elements from the line in line list.                                                                                                                                                                                                                                                       
``method``            str                        ``holy-grail``, ``identify``, ``reidentify``, ``echelle``, ``full_template``  ``holy-grail``    Method to use to fit the individual arc lines.  Note that some of the available methods should not be used; they are unstable and require significant parameter tweaking to succeed.  You should use one of 'holy-grail', 'reidentify', or 'full_template'.  'holy-grail' attempts to get a first guess at line IDs by looking for patterns in the line locations.  It is fully automated.  When it works, it works well; however, it can fail catastrophically.  Instead, 'reidentify' and 'full_template' are the preferred methods.  They require an archived wavelength solution for your specific instrument/grating combination as a reference.  This is used to anchor the wavelength solution for the data being reduced.  All options are: holy-grail, identify, reidentify, echelle, full_template.
``n_final``           int, float, list, ndarray  ..                                                                            4                 Order of final fit to the wavelength solution (there are n_final+1 parameters in the fit). This can be a single number or a list/array providing the value for each slit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``n_first``           int                        ..                                                                            2                 Order of first guess fit to the wavelength solution.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
``nfitpix``           int                        ..                                                                            5                 Number of pixels to fit when deriving the centroid of the arc lines (an odd number is best)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``nlocal_cc``         int                        ..                                                                            11                Size of pixel window used for local cross-correlation computation for each arc line. If not an odd number one will be added to it to make it odd.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
``nreid_min``         int                        ..                                                                            1                 Minimum number of times that a given candidate reidentified line must be properly matched with a line in the arxiv to be considered a good reidentification. If there is a lot of duplication in the arxiv of the spectra in question (i.e. multislit) set this to a number like 1-4. For echelle this depends on the number of solutions in the arxiv.  Set this to 1 for fixed format echelle spectrographs.  For an echelle with a tiltable grating, this will depend on the number of solutions in the arxiv.                                                                                                                                                                                                                                                                                            
``nsnippet``          int                        ..                                                                            2                 Number of spectra to chop the arc spectrum into when ``method`` is 'full_template'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
``numsearch``         int                        ..                                                                            20                Number of brightest arc lines to search for in preliminary identification                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``qa_log``            bool                       ..                                                                            True              Governs whether the wavelength solution arc line QA plots will have log or linear scalingIf True, the scaling will be log, if False linear                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``redo_slits``        int, list                  ..                                                                            ..                Redo the input slit(s) [multislit] or order(s) [echelle]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``reference``         str                        ``arc``, ``sky``, ``pixel``                                                   ``arc``           Perform wavelength calibration with an arc, sky frame.  Use 'pixel' for no wavelength solution.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
``refframe``          str                        ``observed``, ``heliocentric``, ``barycentric``                               ``heliocentric``  Frame of reference for the wavelength calibration.  Options are: observed, heliocentric, barycentric                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
``reid_arxiv``        str                        ..                                                                            ..                Name of the archival wavelength solution file that will be used for the wavelength reidentification.  Only used if ``method`` is 'reidentify' or 'full_template'.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
``rms_threshold``     float                      ..                                                                            0.15              Maximum RMS (in binned pixels) for keeping a slit/order solution. Used for echelle spectrographs, the 'reidentify' method, and when re-analyzing a slit with the redo_slits parameter.In a future PR, we will refactor the code to always scale this threshold off the measured FWHM of the arc lines.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
``sigdetect``         int, float, list, ndarray  ..                                                                            5.0               Sigma threshold above fluctuations for arc-line detection.  Arcs are continuum subtracted and the fluctuations are computed after continuum subtraction.  This can be a single number or a vector (list or numpy array) that provides the detection threshold for each slit.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
``sigrej_final``      float                      ..                                                                            3.0               Number of sigma for rejection for the final guess to the wavelength solution.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
``sigrej_first``      float                      ..                                                                            2.0               Number of sigma for rejection for the first guess to the wavelength solution.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
``use_instr_flag``    bool                       ..                                                                            False             If True, restrict to lines matching the instrument.  WARNING: This is only implemented for shane_kast_red + HolyGrail.  Do not use it unless you really know what you are doing.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``wvrng_arxiv``       list                       ..                                                                            ..                Cut the arxiv template down to this specified wavelength range [min,max]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
====================  =========================  ============================================================================  ================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _coadd1dpar:

Coadd1DPar Keywords
-------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.Coadd1DPar`

====================  ==========  =======  ==========  =========================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                   Type        Options  Default     Description                                                                                                                                                                                                                                                                                                                                                                                                                              
====================  ==========  =======  ==========  =========================================================================================================================================================================================================================================================================================================================================================================================================================================
``chk_version``       bool        ..       True        If True enforce strict PypeIt version checking to ensure that spec1d*.fits files were createdwith the current version of PypeIt                                                                                                                                                                                                                                                                                                          
``coaddfile``         str         ..       ..          Output filename                                                                                                                                                                                                                                                                                                                                                                                                                          
``dloglam``           int, float  ..       ..          Dispersion in units of log10(wave) in case you want to specify it in the get_wave_grid  (for the 'velocity' or 'log10' options), otherwise a median value is computed from the data.                                                                                                                                                                                                                                                     
``dv``                int, float  ..       ..          Dispersion in units of km/s in case you want to specify it in the get_wave_grid  (for the 'velocity' option), otherwise a median value is computed from the data.                                                                                                                                                                                                                                                                        
``dwave``             int, float  ..       ..          Dispersion in Angstroms in case you want to specify it in the get_wave_grid  (for the 'linear' option), otherwise a median value is computed from the data.                                                                                                                                                                                                                                                                              
``ex_value``          str         ..       ``OPT``     The extraction to coadd, i.e. optimal or boxcar. Must be either 'OPT' or 'BOX'                                                                                                                                                                                                                                                                                                                                                           
``filter``            str         ..       ``none``    Filter for scaling.  See flux_calib.load_fitler_file() for naming.  Ignore if none                                                                                                                                                                                                                                                                                                                                                       
``filter_mag``        float       ..       ..          Magnitude of the source in the given filter                                                                                                                                                                                                                                                                                                                                                                                              
``filter_mask``       str, list   ..       ..          List of wavelength regions to mask when doing the scaling (`i.e.`, occasional junk pixels). Colon and comma separateed, e.g.   5552:5559,6010:6030                                                                                                                                                                                                                                                                                       
``flux_value``        bool        ..       True        If True (default), the code will coadd the fluxed spectra (i.e. the FLAM) in the spec1d files. If False, it will coadd the counts.                                                                                                                                                                                                                                                                                                       
``lower``             int, float  ..       3.0         Lower rejection threshold used for rejecting pixels when combining spectra in units of sigma.                                                                                                                                                                                                                                                                                                                                            
``mag_type``          str         ..       ``AB``      Magnitude type.  AB is the only option currently allowed                                                                                                                                                                                                                                                                                                                                                                                 
``maxiter_reject``    int         ..       5           Maximum number of iterations for stacking and rejection. The code stops iterating either when the output mask does not change betweeen successive iterations or when maxiter_reject is reached.                                                                                                                                                                                                                                          
``maxiter_scale``     int         ..       5           Maximum number of iterations performed for rescaling spectra.                                                                                                                                                                                                                                                                                                                                                                            
``maxrej``            int         ..       ..          Coadding performs iterative rejection by comparing each exposure to a preliminary stack of all the exposures. If this parameter is set then it will not reject more than maxrej pixels per iteration of this rejection. The default is None, which means no maximum on rejected pixels.                                                                                                                                                  
``nbests``            list, int   ..       ..          Number of orders to use for estimating the per exposure weights. Default is None, which will just use one fourth of the total number of orders. This is only used for Echelle                                                                                                                                                                                                                                                            
``nmaskedge``         int         ..       2           Number of edge pixels to mask. This should be removed/fixed.                                                                                                                                                                                                                                                                                                                                                                             
``ref_percentile``    int, float  ..       70.0        Percentile used for selecting the minimum SNR cut from a reference spectrum used to robustly determine the median ratio between spectra. This parameter is used by coadd1d.robust_median_ratio as part of the automatic rescaling procedure. Pixels above this percentile cut are deemed the "good" pixels and are used to compute the ratio of two spectra.  This must be a number between 0 and 100.                                   
``scale_method``      str         ..       ``auto``    Method used to rescale the spectra prior to coadding. The options are: 'auto' -- Determine the scaling method automatically based on the S/N ratio which works well.  'poly' -- Polynomial rescaling.  'median' -- Median rescaling  'none' -- Do not rescale.  'hand' -- Pass in hand scaling factors. This option is not well tested.                                                                                                  
``sigrej_exp``        int, float  ..       ..          Rejection threshold used for rejecting exposures with S/N more than sigrej_exp*sigma above the median S/N. If None (the default), no rejection is performed. Currently, only available for multi-slit observations.                                                                                                                                                                                                                      
``sigrej_scale``      int, float  ..       3.0         Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec.                                                                                                                                                                                                                                                                                                                                                    
``sn_clip``           int, float  ..       30.0        Errors are capped during rejection so that the S/N is never greater than sn_clip. This prevents overly aggressive rejection in high S/N ratio spectrum which neverthless differ at a level greater than the formal S/N due to systematics.                                                                                                                                                                                               
``sn_min_medscale``   int, float  ..       0.5         For scale method set to ``auto``, this sets the minimum SNR for which median scaling is attempted.                                                                                                                                                                                                                                                                                                                                       
``sn_min_polyscale``  int, float  ..       2.0         For scale method set to ``auto``, this sets the minimum SNR for which polynomial scaling is attempted.                                                                                                                                                                                                                                                                                                                                   
``sn_smooth_npix``    int, float  ..       ..          Number of pixels to median filter by when computing S/N used to decide how to scale and weight spectra. If set to None (default), the code will determine the effective number of good pixels per spectrum in the stack that is being co-added and use 10% of this neff.                                                                                                                                                                 
``spec_samp_fact``    float       ..       1.0         Make the wavelength grid  sampling finer (spec_samp_fact < 1.0) or coarser (spec_samp_fact > 1.0) by this sampling factor. This basically multiples the 'native' spectral pixels by spec_samp_fact, i.e. units spec_samp_fact are pixels.                                                                                                                                                                                                
``upper``             int, float  ..       3.0         Upper rejection threshold used for rejecting pixels when combining spectra in units of sigma.                                                                                                                                                                                                                                                                                                                                            
``wave_grid_max``     int, float  ..       ..          Used in case you want to specify the maximum wavelength in your wavelength grid, default=None computes from data                                                                                                                                                                                                                                                                                                                         
``wave_grid_min``     int, float  ..       ..          Used in case you want to specify the minimum wavelength in your wavelength grid, default=None computes from data                                                                                                                                                                                                                                                                                                                         
``wave_method``       str         ..       ``linear``  Method used to construct wavelength grid for coadding spectra. The routine that creates the wavelength is :func:`~pypeit.core.wavecal.wvutils.get_wave_grid`. The options are: 'iref' -- Use the first wavelength array.  'velocity' -- Grid is uniform in velocity.  'log10' -- Grid is uniform in log10(wave). This is the same as velocity.  'linear' -- Grid is uniform in lambda.  'concatenate' -- Meld the input wavelength arrays
====================  ==========  =======  ==========  =========================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _coadd2dpar:

Coadd2DPar Keywords
-------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.Coadd2DPar`

====================  =========  =======  ========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                   Type       Options  Default   Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
====================  =========  =======  ========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``exclude_slits``     str, list  ..       ..        Exclude one or more slits from the coaddition. Example syntax -- DET01:175,DET02:205 or MSC02:2234. This and ``only_slits`` are mutually exclusive. If both are provided, ``only_slits`` takes precedence.                                                                                                                                                                                                                                                                                                                                                                                                                  
``manual``            str        ..       ..        Manual extraction parameters. det:spat:spec:fwhm:boxcar_radius. Multiple manual extractions are semi-colon separated, and spat,spec are in the pseudo-image generated by COADD2D.boxcar_radius is optional and in pixels (not arcsec!).                                                                                                                                                                                                                                                                                                                                                                                     
``offsets``           str, list  ..       ``auto``  Offsets for the images being combined (spat pixels). Options are: ``maskdef_offsets``, ``header``, ``auto``, and a list of offsets. Use ``maskdef_offsets`` to use the offsets computed during the slitmask design matching (currently available for these :ref:`slitmask_info_instruments` only). If equal to ``header``, the dither offsets recorded in the header, when available, will be used. If ``auto`` is chosen, PypeIt will try to compute the offsets using a reference object with the highest S/N, or an object selected by the user (see ``user_obj``). If a list of offsets is provided, PypeIt will use it.
``only_slits``        str, list  ..       ..        Restrict coaddition to one or more of slits. Example syntax -- DET01:175,DET02:205 or MSC02:2234. This and ``exclude_slits`` are mutually exclusive. If both are provided, ``only_slits`` takes precedence.                                                                                                                                                                                                                                                                                                                                                                                                                 
``spat_toler``        int        ..       5         This parameter provides the desired tolerance in spatial pixel used to identify slits in different exposures                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
``use_slits4wvgrid``  bool       ..       False     If True, use the slits to set the trace down the center                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``user_obj``          int, list  ..       ..        Object that the user wants to use to compute the weights and/or the offsets for coadding images. For longslit/multislit spectroscopy, provide the ``SLITID`` and the ``OBJID``, separated by comma, of the selected object. For echelle spectroscopy, provide the ``ECH_OBJID`` of the selected object. See :doc:`out_spec1D` for more info about ``SLITID``, ``OBJID`` and ``ECH_OBJID``. If this parameter is not ``None``, it will be used to compute the offsets only if ``offsets = auto``, and it will used to compute the weights only if ``weights = auto``.                                                        
``wave_method``       str        ..       ..        Argument to :func:`~pypeit.core.wavecal.wvutils.get_wave_grid` method, which determines how the 2d coadd wavelength grid is constructed. The default is None, which will use a linear gridfor longslit/multislit coadds and a log10 grid for echelle coadds. Currently supported options with 2d coadding are:* 'iref' -- Use one of the exposures (the first) as the reference for the wavelength grid * 'velocity' -- Grid is uniform in velocity* 'log10'  -- Grid is uniform in log10(wave). This is the same as velocity.* 'linear' -- Grid is uniform in wavelength                                                   
``weights``           str, list  ..       ``auto``  Mode for the weights used to coadd images. Options are: ``auto``, ``uniform``, or a list of weights. If ``auto`` is used, PypeIt will try to compute the weights using a reference object with the highest S/N, or an object selected by the user (see ``user_obj``), if ``uniform`` is used, uniform weights will be applied. If a list of weights is provided, PypeIt will use it.                                                                                                                                                                                                                                        
====================  =========  =======  ========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _collate1dpar:

Collate1DPar Keywords
---------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.Collate1DPar`

=========================  ==========  =======  ============================================  ==================================================================================================================================================================================================================================================================================================================================================================================================================
Key                        Type        Options  Default                                       Description                                                                                                                                                                                                                                                                                                                                                                                                       
=========================  ==========  =======  ============================================  ==================================================================================================================================================================================================================================================================================================================================================================================================================
``chk_version``            bool        ..       False                                         Whether to check the data model versions of spec1d files and sensfunc files.                                                                                                                                                                                                                                                                                                                                      
``dry_run``                bool        ..       False                                         If set, the script will display the matching File and Object Ids but will not flux, coadd or archive.                                                                                                                                                                                                                                                                                                             
``exclude_serendip``       bool        ..       False                                         Whether to exclude SERENDIP objects from collating.                                                                                                                                                                                                                                                                                                                                                               
``exclude_slit_trace_bm``  list, str   ..       []                                            A list of slit trace bitmask bits that should be excluded.                                                                                                                                                                                                                                                                                                                                                        
``flux``                   bool        ..       False                                         If set, the script will flux calibrate using archived sensfuncs before coadding.                                                                                                                                                                                                                                                                                                                                  
``ignore_flux``            bool        ..       False                                         If set, the script will only coadd non-fluxed spectra even if flux data is present. Otherwise fluxed spectra are coadded if all spec1ds have been fluxed calibrated.                                                                                                                                                                                                                                              
``match_using``            str         ..       ``ra/dec``                                    Determines how 1D spectra are matched as being the same object. Must be either 'pixel' or 'ra/dec'.                                                                                                                                                                                                                                                                                                               
``outdir``                 str         ..       ``/Users/westfall/Work/packages/pypeit/doc``  The path where all coadded output files and report files will be placed.                                                                                                                                                                                                                                                                                                                                          
``refframe``               str         ..       ..                                            Perform reference frame correction prior to coadding. Options are: observed, heliocentric, barycentric                                                                                                                                                                                                                                                                                                            
``spec1d_outdir``          str         ..       ..                                            The path where all modified spec1d files are placed. These are only created if flux calibration or refframe correction are asked for.                                                                                                                                                                                                                                                                             
``tolerance``              str, float  ..       ``1.0``                                       The tolerance used when comparing the coordinates of objects. If two objects are within this distance from each other, they are considered the same object. If match_using is 'ra/dec' (the default) this is an angular distance. The defaults units are arcseconds but other units supported by astropy.coordinates.Angle can be used (`e.g.`, '0.003d' or '0h1m30s'). If match_using is 'pixel' this is a float.
``wv_rms_thresh``          float       ..       ..                                            If set, any objects with a wavelength RMS > this value are skipped, else all wavelength RMS values are accepted.                                                                                                                                                                                                                                                                                                  
=========================  ==========  =======  ============================================  ==================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _flexurepar:

FlexurePar Keywords
-------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.FlexurePar`

===================  ==========  ========================================================  ====================  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                  Type        Options                                                   Default               Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
===================  ==========  ========================================================  ====================  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``excessive_shift``  str         ``crash``, ``set_to_zero``, ``continue``, ``use_median``  ``use_median``        Behavior when the measured spectral flexure shift is larger than ``spec_maxshift``.  The options are: 'crash' - Raise an error and halt the data reduction; 'set_to_zero' - Set the flexure shift to zero and continue with the reduction; 'continue' - Use the large flexure value whilst issuing a warning; and 'use_median' - Use the median flexure shift among all the objects in the same slit (if more than one object is detected) or among all the other slits; if not available, the flexure correction will not be applied.
``multi_min_SN``     int, float  ..                                                        1                     Minimum S/N for analyzing sky spectrum for flexure                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``spec_maxshift``    int         ..                                                        20                    Maximum allowed spectral flexure shift in pixels.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``spec_method``      str         ``boxcar``, ``slitcen``, ``skip``                         ``skip``              Method used to correct for flexure. Use skip for no correction.  If slitcen is used, the flexure correction is performed before the extraction of objects (not recommended).  Options are: None, boxcar, slitcen, skip                                                                                                                                                                                                                                                                                                                
``spectrum``         str         ..                                                        ``paranal_sky.fits``  Archive sky spectrum to be used for the flexure correction.                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
===================  ==========  ========================================================  ====================  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _fluxcalibratepar:

FluxCalibratePar Keywords
-------------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.FluxCalibratePar`

=====================  ====  =======  ===========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                    Type  Options  Default      Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
=====================  ====  =======  ===========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``extinct_correct``    bool  ..       ..           The default behavior for atmospheric extinction corrections is that if UVIS algorithm is used (which does not correct for telluric absorption) than an atmospheric extinction model is used to correct for extinction below 10,000A, whereas if the IR algorithm is used, then no extinction correction is applied since the atmosphere is modeled directly. To follow these defaults based on the algorithm this parameter should be set to ``extinct_correct=None``. If instead this parameter is set, this overide this default behavior. In other words, it will force an extinction correction if ``extinct_correct=True``, and will not perform an extinction correction if ``extinct_correct=False``.
``extinct_file``       str   ..       ``closest``  If ``extinct_file='closest'`` the code will select the PypeIt-included extinction file for the closest observatory (within 5 deg, geographic coordinates) to the telescope identified in ``std_file`` (see :ref:`extinction_correction` for the list of currently included files).  If constructing a sesitivity function for a telescope not within 5 deg of a listed observatory, this parameter may be set to the name of one of the listed extinction files.  Alternatively, a custom extinction file may be installed in the PypeIt cache using the ``pypeit_install_extinctfile`` script; this parameter may then be set to the name of the custom extinction file.                                   
``extrap_sens``        bool  ..       False        If False (default), the code will crash if one tries to use sensfunc at wavelengths outside its defined domain. By changing the par['sensfunc']['extrap_blu'] and par['sensfunc']['extrap_red'] this domain can be extended. If True the code will blindly extrapolate.                                                                                                                                                                                                                                                                                                                                                                                                                                     
``use_archived_sens``  bool  ..       False        Use an archived sensfunc to flux calibration                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
=====================  ====  =======  ===========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _reduxpar:

ReduxPar Keywords
-----------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.ReduxPar`

======================  ==============  =======  ============================================  ===============================================================================================================================================================================================================================================================================================================================================================
Key                     Type            Options  Default                                       Description                                                                                                                                                                                                                                                                                                                                                    
======================  ==============  =======  ============================================  ===============================================================================================================================================================================================================================================================================================================================================================
``calwin``              int, float      ..       0                                             The window of time in hours to search for calibration frames for a science frame                                                                                                                                                                                                                                                                               
``detnum``              int, list       ..       ..                                            Restrict reduction to a list of detector indices. In case of mosaic reduction (currently only available for Gemini/GMOS and Keck/DEIMOS) ``detnum`` should be a list of tuples of the detector indices that are mosaiced together. E.g., for Gemini/GMOS ``detnum`` would be ``[(1,2,3)]`` and for Keck/DEIMOS it would be ``[(1, 5), (2, 6), (3, 7), (4, 8)]``
``ignore_bad_headers``  bool            ..       False                                         Ignore bad headers (NOT recommended unless you know it is safe).                                                                                                                                                                                                                                                                                               
``maskIDs``             str, int, list  ..       ..                                            Restrict reduction to a set of slitmask IDs Example syntax -- ``maskIDs = 818006,818015`` This must be used with detnum (for now).                                                                                                                                                                                                                             
``qadir``               str             ..       ``QA``                                        Directory relative to calling directory to write quality assessment files.                                                                                                                                                                                                                                                                                     
``quicklook``           bool            ..       False                                         Run a quick look reduction? This is usually good if you want to quickly reduce the data (usually at the telescope in real time) to get an initial estimate of the data quality.                                                                                                                                                                                
``redux_path``          str             ..       ``/Users/westfall/Work/packages/pypeit/doc``  Path to folder for performing reductions.  Default is the current working directory.                                                                                                                                                                                                                                                                           
``scidir``              str             ..       ``Science``                                   Directory relative to calling directory to write science files.                                                                                                                                                                                                                                                                                                
``slitspatnum``         str, list       ..       ..                                            Restrict reduction to a set of slit DET:SPAT values (closest slit is used). Example syntax -- slitspatnum = DET01:175,DET01:205 or MSC02:2234  If you are re-running the code, (i.e. modifying one slit) you *must* have the precise SPAT_ID index.                                                                                                            
``sortroot``            str             ..       ..                                            A filename given to output the details of the sorted files.  If None, the default is the root name of the pypeit file.  If off, no output is produced.                                                                                                                                                                                                         
``spectrograph``        str             ..       ..                                            Spectrograph that provided the data to be reduced.  See :ref:`instruments` for valid options.                                                                                                                                                                                                                                                                  
======================  ==============  =======  ============================================  ===============================================================================================================================================================================================================================================================================================================================================================


----

.. _reducepar:

ReducePar Keywords
------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.ReducePar`

==============  ============================================  =======  =========================  =================================================================================
Key             Type                                          Options  Default                    Description                                                                      
==============  ============================================  =======  =========================  =================================================================================
``cube``        :class:`~pypeit.par.pypeitpar.CubePar`        ..       `CubePar Keywords`_        Parameters for cube generation algorithms                                        
``extraction``  :class:`~pypeit.par.pypeitpar.ExtractionPar`  ..       `ExtractionPar Keywords`_  Parameters for extraction algorithms                                             
``findobj``     :class:`~pypeit.par.pypeitpar.FindObjPar`     ..       `FindObjPar Keywords`_     Parameters for the find object and tracing algorithms                            
``skysub``      :class:`~pypeit.par.pypeitpar.SkySubPar`      ..       `SkySubPar Keywords`_      Parameters for sky subtraction algorithms                                        
``slitmask``    :class:`~pypeit.par.pypeitpar.SlitMaskPar`    ..       `SlitMaskPar Keywords`_    Parameters for slitmask                                                          
``trim_edge``   list                                          ..       3, 3                       Trim the slit by this number of pixels left/right when performing sky subtraction
==============  ============================================  =======  =========================  =================================================================================


----

.. _cubepar:

CubePar Keywords
----------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.CubePar`

====================  =====  =======  ============  ===========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                   Type   Options  Default       Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
====================  =====  =======  ============  ===========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``align``             bool   ..       False         If set to True, the input frames will be spatially aligned by cross-correlating the whitelight images with either a reference image (see ``reference_image``) or the whitelight image that is generated using the first spec2d listed in the coadd3d file. Alternatively, the user can specify the offsets (i.e. Delta RA x cos(dec) and Delta Dec, both in arcsec) in the spec2d block of the coadd3d file. See the documentation for examples of this usage.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``astrometric``       bool   ..       True          If true, an astrometric correction will be applied using the alignment frames.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``combine``           bool   ..       False         If set to True, the input frames will be combined. Otherwise, a separate datacube will be generated for each input spec2d file, and will be saved as a spec3d file.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``dec_max``           float  ..       ..            Maximum DEC to use when generating the WCS. If None, the default is maximum DEC based on the WCS of all spaxels. Units should be degrees.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``dec_min``           float  ..       ..            Minimum DEC to use when generating the WCS. If None, the default is minimum DEC based on the WCS of all spaxels. Units should be degrees.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``grating_corr``      bool   ..       True          This option performs a small correction for the relative blaze function of all input frames that have (even slightly) different grating angles, or if you are flux calibrating your science data with a standard star that was observed with a slightly different setup.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``method``            str    ..       ``subpixel``  What method should be used to generate the datacube. There are currently two options: (1) "subpixel" (default) - this algorithm divides each pixel in the spec2d frames into subpixels, and assigns each subpixel to a voxel of the datacube. Flux is conserved, but voxels are correlated, and the error spectrum does not account for covariance between adjacent voxels. See also, spec_subpixel and spat_subpixel. (2) "NGP" (nearest grid point) - this algorithm is effectively a 3D histogram. Flux is conserved, voxels are not correlated, however this option suffers the same downsides as any histogram; the choice of bin sizes can change how the datacube appears. This algorithm takes each pixel on the spec2d frame and puts the flux of this pixel into one voxel in the datacube. Depending on the binning used, some voxels may be empty (zero flux) while a neighboring voxel might contain the flux from two spec2d pixels. Note that all spec2d pixels that contribute to the same voxel are inverse variance weighted (e.g. if two pixels have the same variance, the voxel would be assigned the average flux of the two pixels).
``output_filename``   str    ..       ..            If combining multiple frames, this string sets the output filename of the combined datacube. If combine=False, the output filenames will be prefixed with ``spec3d_*``                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``ra_max``            float  ..       ..            Maximum RA to use when generating the WCS. If None, the default is maximum RA based on the WCS of all spaxels. Units should be degrees.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``ra_min``            float  ..       ..            Minimum RA to use when generating the WCS. If None, the default is minimum RA based on the WCS of all spaxels. Units should be degrees.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``reference_image``   str    ..       ..            White light image of a previously combined datacube. The white light image will be used as a reference when calculating the offsets of the input spec2d files. Ideally, the reference image should have the same shape as the data to be combined (i.e. set the ra_min, ra_max etc. params so they are identical to the reference image).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``relative_weights``  bool   ..       False         If set to True, the combined frames will use a relative weighting scheme. This only works well if there is a common continuum source in the field of view of all input observations, and is generally only required if high relative precision is desired.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
``save_whitelight``   bool   ..       False         Save a white light image of the combined datacube. The output filename will be given by the "output_filename" variable with a suffix "_whitelight". Note that the white light image collapses the flux along the wavelength axis, so some spaxels in the 2D white light image may have different wavelength ranges. To set the wavelength range, use the "whitelight_range" parameter. If combine=False, the individual spec3d files will have a suffix "_whitelight".                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``scale_corr``        str    ..       ..            This option performs a small correction for the relative spectral illumination scale of different spec2D files. Specify the relative path+file to the spec2D file that you would like to use for the relative scaling. If you want to perform this correction, it is best to use the spec2d file with the highest S/N sky spectrum. You should choose the same frame for both the standards and science frames.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
``skysub_frame``      str    ..       ``image``     Set the sky subtraction to be implemented. The default behaviour is to subtract the sky using the model that is derived from each individual image (i.e. set this parameter to "image"). To turn off sky subtraction completely, set this parameter to "none" (all lowercase). Finally, if you want to use a different frame for the sky subtraction, specify the relative path+file to the spec2D file that you would like to use for the sky subtraction. The model fit to the sky of the specified frame will be used. Note, the sky and science frames do not need to have the same exposure time; the sky model will be scaled to the science frame based on the relative exposure time.                                                                                                                                                                                                                                                                                                                                                                                                                                                              
``slit_spec``         bool   ..       True          If the data use slits in one spatial direction, set this to True. If the data uses fibres for all spaxels, set this to False.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
``spat_subpixel``     int    ..       5             When method=subpixel, spat_subpixel sets the subpixellation scale of each detector pixel in the spatial direction. The total number of subpixels in each pixel is given by spec_subpixel x spat_subpixel. The default option is to divide each spec2d pixel into 25 subpixels during datacube creation. See also, spec_subpixel.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
``spatial_delta``     float  ..       ..            The spatial size of each spaxel to use when generating the WCS (in arcsec). If None, the default is set by the spectrograph file.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
``spec_subpixel``     int    ..       5             When method=subpixel, spec_subpixel sets the subpixellation scale of each detector pixel in the spectral direction. The total number of subpixels in each pixel is given by spec_subpixel x spat_subpixel. The default option is to divide each spec2d pixel into 25 subpixels during datacube creation. See also, spat_subpixel.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
``standard_cube``     str    ..       ..            Filename of a standard star datacube. This cube will be used to correct the relative scales of the slits, and to flux calibrate the science datacube.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
``wave_delta``        float  ..       ..            The wavelength step to use when generating the WCS (in Angstroms). If None, the default is set by the wavelength solution.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
``wave_max``          float  ..       ..            Maximum wavelength to use when generating the WCS. If None, the default is maximum wavelength based on the WCS of all spaxels. Units should be Angstroms.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``wave_min``          float  ..       ..            Minimum wavelength to use when generating the WCS. If None, the default is minimum wavelength based on the WCS of all spaxels. Units should be Angstroms.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``whitelight_range``  list   ..       None, None    A two element list specifying the wavelength range over which to generate the white light image. The first (second) element is the minimum (maximum) wavelength to use. If either of these elements are None, PypeIt will automatically use a wavelength range that ensures all spaxels have the same wavelength coverage. Note, if you are using a reference_image to align all frames, it is preferable to use the same white light wavelength range for all white light images. For example, you may wish to use an emission line map to register two frames.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
====================  =====  =======  ============  ===========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _extractionpar:

ExtractionPar Keywords
----------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.ExtractionPar`

====================  ==========  =======  =======  =============================================================================================================================================================================================================================================================================================
Key                   Type        Options  Default  Description                                                                                                                                                                                                                                                                                  
====================  ==========  =======  =======  =============================================================================================================================================================================================================================================================================================
``boxcar_radius``     int, float  ..       1.5      Boxcar radius in arcseconds used for boxcar extraction                                                                                                                                                                                                                                       
``model_full_slit``   bool        ..       False    If True local sky subtraction will be performed on the entire slit. If False, local sky subtraction will be applied to only a restricted region around each object. This should be set to True for either multislit observations using narrow slits or echelle observations with narrow slits
``return_negative``   bool        ..       False    If ``True`` the negative traces will be extracted and saved to disk                                                                                                                                                                                                                          
``skip_extraction``   bool        ..       False    Do not perform an object extraction                                                                                                                                                                                                                                                          
``skip_optimal``      bool        ..       False    Perform boxcar extraction only (i.e. skip Optimal and local skysub)                                                                                                                                                                                                                          
``sn_gauss``          int, float  ..       4.0      S/N threshold for performing the more sophisticated optimal extraction which performs a b-spline fit to the object profile. For S/N < sn_gauss the code will simply optimal extractwith a Gaussian with FWHM determined from the object finding.                                             
``std_prof_nsigma``   float       ..       30.0     prof_nsigma parameter for Standard star extraction.  Prevents undesired rejection. NOTE: Not consumed by the code at present.                                                                                                                                                                
``use_2dmodel_mask``  bool        ..       True     Mask pixels rejected during profile fitting when extracting.Turning this off may help with bright emission lines.                                                                                                                                                                            
``use_user_fwhm``     bool        ..       False    Boolean indicating if PypeIt should use the FWHM provided by the user (``find_fwhm`` in `FindObjPar`) for the optimal extraction. If this parameter is ``False`` (default), PypeIt estimates the FWHM for each detected object, and uses ``find_fwhm`` as initial guess.                     
====================  ==========  =======  =======  =============================================================================================================================================================================================================================================================================================


----

.. _findobjpar:

FindObjPar Keywords
-------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.FindObjPar`

===========================  ==========  =======  =======  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                          Type        Options  Default  Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
===========================  ==========  =======  =======  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``ech_find_max_snr``         int, float  ..       1.0      Criteria for keeping echelle objects. They must either have a maximum S/N across all the orders greater than this value  or satisfy the min_snr criteria described by the min_snr parameters. If maxnumber is set (see above) then these criteria will be applied but only the maxnumber highest (median) S/N ratio objects will be kept.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``ech_find_min_snr``         int, float  ..       0.3      Criteria for keeping echelle objects. They must either have a maximum S/N across all the orders greater than ech_find_max_snr,  value  or they must have S/N > ech_find_min_snr on >= ech_find_nabove_min_snr orders. If maxnumber is set (see above) then these criteria will be applied but only the maxnumber highest (median) S/N ratio objects will be kept.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
``ech_find_nabove_min_snr``  int         ..       2        Criteria for keeping echelle objects. They must either have a maximum S/N across all the orders greater than ech_find_max_snr,  value  or they must have S/N > ech_find_min_snr on >= ech_find_nabove_min_snr orders. If maxnumber is set (see above) then these criteria will be applied but only the maxnumber highest (median) S/N ratio objects will be kept.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
``find_extrap_npoly``        int         ..       3        Polynomial order used for trace extrapolation                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
``find_fwhm``                int, float  ..       5.0      Indicates roughly the fwhm of objects in pixels for object finding                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
``find_maxdev``              int, float  ..       2.0      Maximum deviation of pixels from polynomial fit to trace used to reject bad pixels in trace fitting.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
``find_min_max``             list        ..       ..       It defines the minimum and maximum of your object in pixels in the spectral direction on the detector. It only used for object finding. This parameter is helpful if your object only has emission lines or at high redshift and the trace only shows in part of the detector.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
``find_negative``            bool        ..       ..       Identify negative objects in object finding for spectra that are differenced. This is used to manually override the default behavior in PypeIt for object finding by setting this parameter to something other than None The default behavior is that PypeIt will search for negative object traces if background frames are present in the PypeIt file that are classified as "science" (i.e. via pypeit_setup -b, and setting bkg_id in the PypeIt file). If background frames are present that are classified as "sky", then PypeIt will NOT search for negative object traces. If one wishes to explicitly override this default behavior, set this parameter to True to find negative objects or False to ignore them.                                                                                                                                                                                                                                                                                                                                                                                                                                  
``find_trim_edge``           list        ..       5, 5     Trim the slit by this number of pixels left/right before finding objects                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``maxnumber_sci``            int         ..       10       Maximum number of objects to extract in a science frame.  Use None for no limit. This parameter can be useful in situations where systematics lead to spurious extra objects. Setting this parameter means they will be trimmed. For mulitslit maxnumber applies per slit, for echelle observations this applies per order. Note that objects on a slit/order impact the sky-modeling and so maxnumber should never be lower than the true number of detectable objects on your slit. For image differenced observations with positive and negative object traces, maxnumber applies to the number of positive (or negative) traces individually. In other words, if you had two positive objects and one negative object, then you would set maxnumber to be equal to two (not three). Note that if manually extracted apertures are explicitly requested, they do not count against this maxnumber. If more than maxnumber objects are detected, then highest S/N ratio objects will be the ones that are kept. For multislit observations the choice here depends on the slit length. For echelle observations with short slits we set the default to be 1
``maxnumber_std``            int         ..       5        Maximum number of objects to extract in a standard star frame.  Same functionality as maxnumber_sci documented above. For multislit observations the default here is 5, for echelle observations the default is 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
``skip_final_global``        bool        ..       False    If True, do not update initial sky to get global sky using updated noise model. This should be True for quicklook to save time. This should also be True for near-IR reductions which perform difference imaging, since there we fit sky-residuals rather than the sky itself, so there is no noise model to update.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
``skip_second_find``         bool        ..       False    Only perform one round of object finding (mainly for quick_look)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``skip_skysub``              bool        ..       False    If True, do not sky subtract when performing object finding. This should be set to True for example when running on data that is already sky-subtracted. Note that for near-IR difference imaging one still wants to remove sky-residuals via sky-subtraction, and so this is typically set to False                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
``snr_thresh``               int, float  ..       10.0     S/N threshold for object finding in wavelength direction smashed image.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
``std_spec1d``               str         ..       ..       A PypeIt spec1d file of a previously reduced standard star.  The trace of the standard star spectrum is used as a crutch for tracing the object spectra, when a direct trace is not possible (i.e., faint sources).  If provided, this overrides use of any standards included in your pypeit file; the standard exposures will still be reduced.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
``trace_npoly``              int         ..       5        Order of legendre polynomial fits to object traces.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
===========================  ==========  =======  =======  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _skysubpar:

SkySubPar Keywords
------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.SkySubPar`

===================  ==========  =======  =======  ===================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                  Type        Options  Default  Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
===================  ==========  =======  =======  ===================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``bspline_spacing``  int, float  ..       0.6      Break-point spacing for the bspline sky subtraction fits.                                                                                                                                                                                                                                                                                                                                                                                                                          
``global_sky_std``   bool        ..       True     Global sky subtraction will be performed on standard stars. This should be turned off for example for near-IR reductions with narrow slits, since bright standards can fill the slit causing global sky-subtraction to fail. In these situations we go straight to local sky-subtraction since it is designed to deal with such situations                                                                                                                                         
``joint_fit``        bool        ..       False    Perform a simultaneous joint fit to sky regions using all available slits. Currently, this parameter is only used for IFU data reduction. Note that the current implementation does not account for variations in the instrument FWHM in different slits. This will be addressed by Issue #1660.                                                                                                                                                                                   
``local_maskwidth``  float       ..       4.0      Initial width of the region in units of FWHM that will be used for local sky subtraction                                                                                                                                                                                                                                                                                                                                                                                           
``mask_by_boxcar``   bool        ..       False    In global sky evaluation, mask the sky region around the object by the boxcar radius (set in ExtractionPar).                                                                                                                                                                                                                                                                                                                                                                       
``max_mask_frac``    float       ..       0.8      Maximum fraction of total pixels on a slit that can be masked by the input masks. If more than this threshold is masked the code will return zeros and throw a warning.                                                                                                                                                                                                                                                                                                            
``no_local_sky``     bool        ..       False    If True, turn off local sky model evaluation, but do fit object profile and perform optimal extraction                                                                                                                                                                                                                                                                                                                                                                             
``no_poly``          bool        ..       False    Turn off polynomial basis (Legendre) in global sky subtraction                                                                                                                                                                                                                                                                                                                                                                                                                     
``sky_sigrej``       float       ..       3.0      Rejection parameter for local sky subtraction                                                                                                                                                                                                                                                                                                                                                                                                                                      
``user_regions``     str, list   ..       ..       Provides a user-defined mask defining sky regions.  By default, the sky regions are identified automatically.  To specify sky regions for *all* slits, provide a comma separated list of percentages.  For example, setting user_regions = :10,35:65,80: selects the first 10%, the inner 30%, and the final 20% of *all* slits as containing sky.  Setting user_regions = user will attempt to load any SkyRegions files generated by the user via the pypeit_skysub_regions tool.
===================  ==========  =======  =======  ===================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _slitmaskpar:

SlitMaskPar Keywords
--------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.SlitMaskPar`

===========================  ==========  =======  =======  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                          Type        Options  Default  Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
===========================  ==========  =======  =======  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``assign_obj``               bool        ..       False    If SlitMask object was generated, assign RA,DEC,name to detected objects                                                                                                                                                                                                                                                                                                                                                                                                                                              
``bright_maskdef_id``        int         ..       ..       `maskdef_id` (corresponding e.g., to `dSlitId` and `Slit_Number` in the DEIMOS/LRIS and MOSFIRE slitmask design, respectively) of a slit containing a bright object that will be used to compute the slitmask offset. This parameter is optional and is ignored if ``slitmask_offset`` is provided.                                                                                                                                                                                                                   
``extract_missing_objs``     bool        ..       False    Force extraction of undetected objects at the location expected from the slitmask design.                                                                                                                                                                                                                                                                                                                                                                                                                             
``missing_objs_boxcar_rad``  int, float  ..       1.0      Indicates the boxcar radius in arcsec for the force extraction of undetected objects.                                                                                                                                                                                                                                                                                                                                                                                                                                 
``missing_objs_fwhm``        int, float  ..       ..       Indicates the FWHM in arcsec for the force extraction of undetected objects. PypeIt will try to determine the FWHM from the flux profile (by using ``missing_objs_fwhm`` as initial guess). If the FWHM cannot be determined, ``missing_objs_fwhm`` will be assumed. If you do not want PypeIt to try to determine the FWHM set the parameter ``use_user_fwhm`` in ``ExtractionPar`` to True. If ``missing_objs_fwhm`` is ``None`` (which is the default) PypeIt will use the median FWHM of all the detected objects.
``obj_toler``                int, float  ..       1.0      If slitmask design information is provided, and slit matching is performed (``use_maskdesign = True`` in ``EdgeTracePar``), this parameter provides the desired tolerance (arcsec) to match sources to targeted objects                                                                                                                                                                                                                                                                                               
``slitmask_offset``          int, float  ..       ..       User-provided slitmask offset (pixels) from the position expected by the slitmask design. This is optional, and if set PypeIt will NOT compute the offset using ``snr_thrshd`` or ``bright_maskdef_id``.                                                                                                                                                                                                                                                                                                              
``snr_thrshd``               int, float  ..       50.0     Objects detected above this S/N threshold will be used to compute the slitmask offset. This is the default behaviour for DEIMOS  unless ``slitmask_offset``, ``bright_maskdef_id`` or ``use_alignbox`` is set.                                                                                                                                                                                                                                                                                                        
``use_alignbox``             bool        ..       False    Use stars in alignment boxes to compute the slitmask offset. If this is set to ``True`` PypeIt will NOT compute the offset using ``snr_thrshd`` or ``bright_maskdef_id``                                                                                                                                                                                                                                                                                                                                              
``use_dither_offset``        bool        ..       False    Use the dither offset recorded in the header of science frames as the value of the slitmask offset. This is currently only available for Keck MOSFIRE reduction and it is set as the default for this instrument. If set PypeIt will NOT compute the offset using ``snr_thrshd`` or ``bright_maskdef_id``. However, it is ignored if ``slitmask_offset`` is provided.                                                                                                                                                 
===========================  ==========  =======  =======  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _framegrouppar:

FrameGroupPar Keywords
----------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.FrameGroupPar`

=============  ===============================================  ============================================================================================================================================================  ============================  ===============================================================================================================================================================================================================================================================
Key            Type                                             Options                                                                                                                                                       Default                       Description                                                                                                                                                                                                                                                    
=============  ===============================================  ============================================================================================================================================================  ============================  ===============================================================================================================================================================================================================================================================
``exprng``     list                                             ..                                                                                                                                                            None, None                    Used in identifying frames of this type.  This sets the minimum and maximum allowed exposure times.  There must be two items in the list.  Use None to indicate no limit; i.e., to select exposures with any time greater than 30 sec, use exprng = [30, None].
``frametype``  str                                              ``align``, ``arc``, ``bias``, ``dark``, ``pinhole``, ``pixelflat``, ``illumflat``, ``lampoffflats``, ``science``, ``standard``, ``trace``, ``tilt``, ``sky``  ``science``                   Frame type.  Options are: align, arc, bias, dark, pinhole, pixelflat, illumflat, lampoffflats, science, standard, trace, tilt, sky                                                                                                                             
``process``    :class:`~pypeit.par.pypeitpar.ProcessImagesPar`  ..                                                                                                                                                            `ProcessImagesPar Keywords`_  Low level parameters used for basic image processing                                                                                                                                                                                                           
``useframe``   str                                              ..                                                                                                                                                            ..                            A calibrations file to use if it exists.                                                                                                                                                                                                                       
=============  ===============================================  ============================================================================================================================================================  ============================  ===============================================================================================================================================================================================================================================================


----

.. _processimagespar:

ProcessImagesPar Keywords
-------------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.ProcessImagesPar`

========================  ==========  ======================================  ==========  ============================================================================================================================================================================================================================================================================================================================================================
Key                       Type        Options                                 Default     Description                                                                                                                                                                                                                                                                                                                                                 
========================  ==========  ======================================  ==========  ============================================================================================================================================================================================================================================================================================================================================================
``apply_gain``            bool        ..                                      True        Convert the ADUs to electrons using the detector gain                                                                                                                                                                                                                                                                                                       
``clip``                  bool        ..                                      True        Perform sigma clipping when combining.  Only used with combine=mean                                                                                                                                                                                                                                                                                         
``comb_sigrej``           float       ..                                      ..          Sigma-clipping level for when clip=True; Use None for automatic limit (recommended).                                                                                                                                                                                                                                                                        
``combine``               str         ``median``, ``mean``                    ``mean``    Method used to combine multiple frames.  Options are: median, mean                                                                                                                                                                                                                                                                                          
``dark_expscale``         bool        ..                                      False       If designated dark frames are used and have a different exposure time than the science frames, scale the counts by the by the ratio in the exposure times to adjust the dark counts for the difference in exposure time.  WARNING: You should always take dark frames that have the same exposure time as your science frames, so use this option with care!
``empirical_rn``          bool        ..                                      False       If True, use the standard deviation in the overscan region to measure an empirical readnoise to use in the noise model.                                                                                                                                                                                                                                     
``grow``                  int, float  ..                                      1.5         Factor by which to expand regions with cosmic rays detected by the LA cosmics routine.                                                                                                                                                                                                                                                                      
``lamaxiter``             int         ..                                      1           Maximum number of iterations for LA cosmics routine.                                                                                                                                                                                                                                                                                                        
``mask_cr``               bool        ..                                      False       Identify CRs and mask them                                                                                                                                                                                                                                                                                                                                  
``n_lohi``                list        ..                                      0, 0        Number of pixels to reject at the lowest and highest ends of the distribution; i.e., n_lohi = low, high.  Use None for no limit.                                                                                                                                                                                                                            
``noise_floor``           float       ..                                      0.0         Impose a noise floor by adding the provided fraction of the bias- and dark-subtracted electron counts to the error budget.  E.g., a value of 0.01 means that the S/N of the counts in the image will never be greater than 100.                                                                                                                             
``objlim``                int, float  ..                                      3.0         Object detection limit in LA cosmics routine                                                                                                                                                                                                                                                                                                                
``orient``                bool        ..                                      True        Orient the raw image into the PypeIt frame                                                                                                                                                                                                                                                                                                                  
``overscan_method``       str         ``polynomial``, ``savgol``, ``median``  ``savgol``  Method used to fit the overscan. Options are: polynomial, savgol, median                                                                                                                                                                                                                                                                                    
``overscan_par``          int, list   ..                                      5, 65       Parameters for the overscan subtraction.  For 'polynomial', set overcan_par = order, number of pixels, number of repeats ; for 'savgol', set overscan_par = order, window size ; for 'median', set overscan_par = None or omit the keyword.                                                                                                                 
``rmcompact``             bool        ..                                      True        Remove compact detections in LA cosmics routine                                                                                                                                                                                                                                                                                                             
``satpix``                str         ``reject``, ``force``, ``nothing``      ``reject``  Handling of saturated pixels.  Options are: reject, force, nothing                                                                                                                                                                                                                                                                                          
``shot_noise``            bool        ..                                      True        Use the bias- and dark-subtracted image to calculate and include electron count shot noise in the image processing error budget                                                                                                                                                                                                                             
``sigclip``               int, float  ..                                      4.5         Sigma level for rejection in LA cosmics routine                                                                                                                                                                                                                                                                                                             
``sigfrac``               int, float  ..                                      0.3         Fraction for the lower clipping threshold in LA cosmics routine.                                                                                                                                                                                                                                                                                            
``spat_flexure_correct``  bool        ..                                      False       Correct slits, illumination flat, etc. for flexure                                                                                                                                                                                                                                                                                                          
``subtract_continuum``    bool        ..                                      False       Subtract off the continuum level from an image. This parameter should only be set to True to combine arcs with multiple different lamps. For all other cases, this parameter should probably be False.                                                                                                                                                      
``trim``                  bool        ..                                      True        Trim the image to the detector supplied region                                                                                                                                                                                                                                                                                                              
``use_biasimage``         bool        ..                                      True        Use a bias image.  If True, one or more must be supplied in the PypeIt file.                                                                                                                                                                                                                                                                                
``use_darkimage``         bool        ..                                      False       Subtract off a dark image.  If True, one or more darks must be provided.                                                                                                                                                                                                                                                                                    
``use_illumflat``         bool        ..                                      True        Use the illumination flat to correct for the illumination profile of each slit.                                                                                                                                                                                                                                                                             
``use_overscan``          bool        ..                                      True        Subtract off the overscan.  Detector *must* have one or code will crash.                                                                                                                                                                                                                                                                                    
``use_pattern``           bool        ..                                      False       Subtract off a detector pattern. This pattern is assumed to be sinusoidal along one direction, with a frequency that is constant across the detector.                                                                                                                                                                                                       
``use_pixelflat``         bool        ..                                      True        Use the pixel flat to make pixel-level corrections.  A pixelflat image must be provied.                                                                                                                                                                                                                                                                     
``use_specillum``         bool        ..                                      False       Use the relative spectral illumination profiles to correct the spectral illumination profile of each slit. This is primarily used for IFUs.  To use this, you must set ``slit_illum_relative=True`` in the ``flatfield`` parameter set!                                                                                                                     
========================  ==========  ======================================  ==========  ============================================================================================================================================================================================================================================================================================================================================================


----

.. _sensfuncpar:

SensFuncPar Keywords
--------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.SensFuncPar`

=======================  ==============================================  ================  ===========================  ============================================================================================================================================================================================================================================================================================================================================================================================
Key                      Type                                            Options           Default                      Description                                                                                                                                                                                                                                                                                                                                                                                 
=======================  ==============================================  ================  ===========================  ============================================================================================================================================================================================================================================================================================================================================================================================
``IR``                   :class:`~pypeit.par.pypeitpar.TelluricPar`      ..                `TelluricPar Keywords`_      Parameters for the IR sensfunc algorithm                                                                                                                                                                                                                                                                                                                                                    
``UVIS``                 :class:`~pypeit.par.pypeitpar.SensfuncUVISPar`  ..                `SensfuncUVISPar Keywords`_  Parameters for the UVIS sensfunc algorithm                                                                                                                                                                                                                                                                                                                                                  
``algorithm``            str                                             ``UVIS``, ``IR``  ``UVIS``                     Specify the algorithm for computing the sensitivity function. The options are:  (1) UVIS = Should be used for data with :math:`\lambda < 7000` A. No detailed model of telluric absorption but corrects for atmospheric extinction.  (2) IR = Should be used for data with :math:`\lambda > 7000` A. Peforms joint fit for sensitivity function and telluric absorption using HITRAN models.
``extrap_blu``           float                                           ..                0.1                          Fraction of minimum wavelength coverage to grow the wavelength coverage of the sensitivitity function in the blue direction (`i.e.`, if the standard star spectrum cuts off at ``wave_min``) the sensfunc will be extrapolated to cover down to  (1.0 - ``extrap_blu``) * ``wave_min``                                                                                                      
``extrap_red``           float                                           ..                0.1                          Fraction of maximum wavelength coverage to grow the wavelength coverage of the sensitivitity function in the red direction (`i.e.`, if the standard star spectrumcuts off at ``wave_max``) the sensfunc will be extrapolated to cover up to  (1.0 + ``extrap_red``) * ``wave_max``                                                                                                          
``flatfile``             str                                             ..                ..                           Flat field file to be used if the sensitivity function model will utilize the blaze function computed from a flat field file in the Calibrations directory, e.g.Calibrations/Flat_A_0_DET01.fits                                                                                                                                                                                            
``hydrogen_mask_wid``    float                                           ..                10.0                         Mask width from line center for hydrogen recombination lines in Angstroms (total mask width is 2x this value).                                                                                                                                                                                                                                                                              
``mask_helium_lines``    bool                                            ..                False                        Mask certain ``HeII`` recombination lines prominent in O-type stars in the sensitivity function fit A region equal to 0.5 * ``hydrogen_mask_wid`` on either side of the line center is masked.                                                                                                                                                                                              
``mask_hydrogen_lines``  bool                                            ..                True                         Mask hydrogen Balmer, Paschen, Brackett, and Pfund recombination lines in the sensitivity function fit. A region equal to ``hydrogen_mask_wid`` on either side of the line center is masked.                                                                                                                                                                                                
``multi_spec_det``       list                                            ..                ..                           List of detectors (identified by their string name, like DET01) to splice together for multi-detector instruments (e.g. DEIMOS). It is assumed that there is *no* overlap in wavelength across detectors (might be ok if there is).  If entered as a list of integers, they should be converted to the detector name.  **Cannot be used with detector mosaics.**                            
``polyorder``            int, list                                       ..                5                            Polynomial order for sensitivity function fitting                                                                                                                                                                                                                                                                                                                                           
``samp_fact``            float                                           ..                1.5                          Sampling factor to make the wavelength grid for sensitivity function finer or coarser. samp_fact > 1.0 oversamples (finer), samp_fact < 1.0 undersamples (coarser).                                                                                                                                                                                                                         
``star_dec``             float                                           ..                ..                           DEC of the standard star. This will override values in the header (`i.e.`, if they are wrong or absent)                                                                                                                                                                                                                                                                                     
``star_mag``             float                                           ..                ..                           Magnitude of the standard star (for near-IR mainly)                                                                                                                                                                                                                                                                                                                                         
``star_ra``              float                                           ..                ..                           RA of the standard star. This will override values in the header (`i.e.`, if they are wrong or absent)                                                                                                                                                                                                                                                                                      
``star_type``            str                                             ..                ..                           Spectral type of the standard star (for near-IR mainly)                                                                                                                                                                                                                                                                                                                                     
=======================  ==============================================  ================  ===========================  ============================================================================================================================================================================================================================================================================================================================================================================================


----

.. _sensfuncuvispar:

SensfuncUVISPar Keywords
------------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.SensfuncUVISPar`

====================  ==========  =======  ===========  =========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                   Type        Options  Default      Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
====================  ==========  =======  ===========  =========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``extinct_correct``   bool        ..       True         If ``extinct_correct=True`` the code will use an atmospheric extinction model to extinction correct the data below 10000A. Note that this correction makes no sense if one is telluric correcting and this shold be set to False                                                                                                                                                                                                                                                                                                                                                                                                                                         
``extinct_file``      str         ..       ``closest``  If ``extinct_file='closest'`` the code will select the PypeIt-included extinction file for the closest observatory (within 5 deg, geographic coordinates) to the telescope identified in ``std_file`` (see :ref:`extinction_correction` for the list of currently included files).  If constructing a sesitivity function for a telescope not within 5 deg of a listed observatory, this parameter may be set to the name of one of the listed extinction files.  Alternatively, a custom extinction file may be installed in the PypeIt cache using the ``pypeit_install_extinctfile`` script; this parameter may then be set to the name of the custom extinction file.
``nresln``            int, float  ..       20           Parameter governing the spacing of the bspline breakpoints in terms of number of resolution elements.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``polycorrect``       bool        ..       True         Whether you want to correct the sensfunc with polynomial in the telluric and recombination line regions                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``polyfunc``          bool        ..       False        Whether you want to use the polynomial fit as your final SENSFUNC                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``resolution``        int, float  ..       3000.0       Expected resolution of the standard star spectrum. This should be measured from the data.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
``sensfunc``          str         ..       ..           FITS file that contains or will contain the sensitivity function.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``std_file``          str         ..       ..           Standard star file to generate sensfunc                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``std_obj_id``        str, int    ..       ..           Specifies object in spec1d file to use as standard. The brightest object found is used otherwise.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``telluric``          bool        ..       False        If ``telluric=True`` the code creates a synthetic standard star spectrum using the Kurucz models, the sens func is created setting nresln=1.5 it contains the correction for telluric lines.                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``telluric_correct``  bool        ..       False        If ``telluric_correct=True`` the code will grab the sens_dict['telluric'] tag from the sensfunc dictionary and apply it to the data.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``trans_thresh``      float       ..       0.9          Parameter for selecting telluric regions which are masked. Locations below this transmission value are masked. If you have significant telluric absorption you should be using telluric.sensnfunc_telluric                                                                                                                                                                                                                                                                                                                                                                                                                                                               
====================  ==========  =======  ===========  =========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


----

.. _telluricpar:

TelluricPar Keywords
--------------------

Class Instantiation: :class:`~pypeit.par.pypeitpar.TelluricPar`

=======================  ==================  =======  ==========================  =================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Key                      Type                Options  Default                     Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
=======================  ==================  =======  ==========================  =================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
``bal_wv_min_max``       list, ndarray       ..       ..                          Min/max wavelength of broad absorption features. If there are several BAL features, the format for this mask is [wave_min_bal1, wave_max_bal1,wave_min_bal2, wave_max_bal2,...]. These masked pixels will be ignored during the fitting.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
``bounds_norm``          tuple               ..       (0.1, 3.0)                  Normalization bounds for scaling the initial object model.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
``delta_coeff_bounds``   tuple               ..       (-20.0, 20.0)               Parameters setting the polynomial coefficient bounds for sensfunc optimization.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
``delta_redshift``       float               ..       0.1                         Range within the redshift can be varied for telluric fitting, i.e. the code performs a bounded optimization within the redshift +- delta_redshift                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
``disp``                 bool                ..       False                       Argument for scipy.optimize.differential_evolution which will display status messages to the screen indicating the status of the optimization. See documentation for telluric.Telluric for a description of the output and how to know if things are working well.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
``fit_wv_min_max``       list                ..       ..                          Pixels within this mask will be used during the fitting. The format is the same with bal_wv_min_max, but this mask is good pixel masks.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
``func``                 str                 ..       ``legendre``                Polynomial model function                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``lower``                int, float          ..       3.0                         Lower rejection threshold in units of sigma_corr*sigma, where sigma is the formal noise of the spectrum, and sigma_corr is an empirically determined correction to the formal error. The distribution of input chi (defined by chi = (data - model)/sigma) values is analyzed, and a correction factor to the formal error sigma_corr is returned which is multiplied into the formal errors. In this way, a rejection threshold of i.e. 3-sigma, will always correspond to roughly the same percentile.  This renormalization is performed with coadd1d.renormalize_errors function, and guarantees that rejection is not too agressive in cases where the empirical errors determined from the chi-distribution differ significantly from the formal noise which is used to determine chi.                                                                                                                     
``mask_lyman_a``         bool                ..       True                        Mask the blueward of Lyman-alpha line during the fitting?                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``maxiter``              int                 ..       2                           Maximum number of iterations for the telluric + object model fitting. The code performs multiple iterations rejecting outliers at each step. The fit is then performed anew to the remaining good pixels. For this reason if you run with the disp=True option, you will see that the f(x) loss function gets progressively better during the iterations.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
``minmax_coeff_bounds``  tuple               ..       (-5.0, 5.0)                 Parameters setting the polynomial coefficient bounds for sensfunc optimization. Bounds are currently determined as follows. We compute an initial fit to the sensfunc in the :func:`~pypeit.core.telluric.init_sensfunc_model` function. That deterines a set of coefficients. The bounds are then determined according to: [(np.fmin(np.abs(this_coeff)*obj_params['delta_coeff_bounds'][0], obj_params['minmax_coeff_bounds'][0]), np.fmax(np.abs(this_coeff)*obj_params['delta_coeff_bounds'][1], obj_params['minmax_coeff_bounds'][1]))]                                                                                                                                                                                                                                                                                                                                                                     
``model``                str                 ..       ``exp``                     Types of polynomial model. Options are poly, square, exp corresponding to normal polynomial, squared polynomial, or exponentiated polynomial                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``npca``                 int                 ..       8                           Number of pca for the objmodel=qso qso PCA fit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``objmodel``             str                 ..       ..                          The object model to be used for telluric fitting. Currently the options are: qso, star, and poly                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
``only_orders``          int, list, ndarray  ..       ..                          Order number, or list of order numbers if you only want to fit specific orders                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``pca_file``             str                 ..       ``qso_pca_1200_3100.fits``  Fits file containing quasar PCA model. Needed for objmodel=qso.  NOTE: This parameter no longer includes the full pathname to the Telluric Model file, but is just the filename of the model itself.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``pca_lower``            int, float          ..       1220.0                      Minimum wavelength for the qso pca model                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
``pca_upper``            int, float          ..       3100.0                      Maximum wavelength for the qso pca model                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
``pix_shift_bounds``     tuple               ..       (-5.0, 5.0)                 Bounds for the pixel shift optimization in telluric model fit in units of pixels. The atmosphere will be allowed to shift within this range during the fit.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
``polish``               bool                ..       True                        If True then differential evolution will perform an additional optimizatino at the end to polish the best fit at the end, which can improve the optimization slightly. See scipy.optimize.differential_evolution for details.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``polyorder``            int                 ..       3                           Order of the polynomial model fit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
``popsize``              int                 ..       30                          A multiplier for setting the total population size for the differential evolution optimization. See scipy.optimize.differential_evolution for details.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
``recombination``        int, float          ..       0.7                         The recombination constant for the differential evolution optimization. This should be in the range [0, 1]. See scipy.optimize.differential_evolution for details.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
``redshift``             int, float          ..       0.0                         The redshift for the object model. This is currently only used by objmodel=qso                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
``resln_frac_bounds``    tuple               ..       (0.5, 1.5)                  Bounds for the resolution fit optimization which is part of the telluric model. This range is in units of the resln_guess, so the (0.5, 1.5) would bound the spectral resolution fit to be within the range bounds_resln = (0.5*resln_guess, 1.5*resln_guess)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
``resln_guess``          int, float          ..       ..                          A guess for the resolution of your spectrum expressed as lambda/dlambda. The resolution is fit explicitly as part of the telluric model fitting, but this guess helps determine the bounds for the optimization (see next). If not provided, the  wavelength sampling of your spectrum will be used and the resolution calculated using a typical sampling of 3 spectral pixels per resolution element.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
``seed``                 int                 ..       777                         An initial seed for the differential evolution optimization, which is a random process. The default is a seed = 777 which will be used to generate a unique seed for every order. A specific seed is used because otherwise the random number generator will use the time for the seed, and the results will not be reproducible.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
``sn_clip``              int, float          ..       30.0                        This adds an error floor to the ivar, preventing too much rejection at high-S/N (`i.e.`, standard stars, bright objects) using the function utils.clip_ivar. A small erorr is added to the input ivar so that the output ivar_out will never give S/N greater than sn_clip. This prevents overly aggressive rejection in high S/N ratio spectra which neverthless differ at a level greater than the formal S/N due to the fact that our telluric models are only good to about 3%.                                                                                                                                                                                                                                                                                                                                                                                                                              
``star_dec``             float               ..       ..                          Object declination in decimal deg                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
``star_mag``             float, int          ..       ..                          AB magnitude in V band                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
``star_ra``              float               ..       ..                          Object right-ascension in decimal deg                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
``star_type``            str                 ..       ..                          stellar type                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
``sticky``               bool                ..       True                        Sticky parameter for the utils.djs_reject algorithm for iterative model fit rejection.  If set to True then points rejected from a previous iteration are kept rejected, in other words the bad pixel mask is the OR of all previous iterations and rejected pixels accumulate. If set to False, the bad pixel mask is the mask from the previous iteration, and if the model fit changes between iterations, points can alternate from being rejected to not rejected. At present this code only performs optimizations with differential evolution and experience shows that sticky needs to be True in order for these to converge. This is because the outliers can be so large that they dominate the loss function, and one never iteratively converges to a good model fit. In other words, the deformations in the model between iterations with sticky=False are too small to approach a reasonable fit.
``telgridfile``          str                 ..       ..                          File containing the telluric grid for the observatory in question. These grids are generated from HITRAN models for each observatory using nominal site parameters. They must be downloaded from the GoogleDrive and installed in your PypeIt installation via the pypeit_install_telluric script. NOTE: This parameter no longer includes the full pathname to the Telluric Grid file, but is just the filename of the grid itself.                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
``tell_norm_thresh``     int, float          ..       0.9                         Threshold of telluric absorption region                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
``tol``                  float               ..       0.001                       Relative tolerance for converage of the differential evolution optimization. See scipy.optimize.differential_evolution for details.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
``upper``                int, float          ..       3.0                         Upper rejection threshold in units of sigma_corr*sigma, where sigma is the formal noise of the spectrum, and sigma_corr is an empirically determined correction to the formal error. See above for description.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
=======================  ==================  =======  ==========================  =================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================



.. _instr_par:

Instrument-Specific Default Configuration
=========================================

The following provides the changes to the global default parameters
provided above for each instrument.  That is, if one were to include
these in the PypeIt file, you would be reproducing the effect of the
`default_pypeit_par` method specific to each derived
:class:`~pypeit.spectrographs.spectrograph.Spectrograph` class.

.. _instr_par-bok_bc:

BOK BC (``bok_bc``)
-------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = bok_bc
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = None, 120,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          exprng = None, 120,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[wavelengths]]
          lamps = NeI, ArI, ArII, HeI,
          fwhm = 5.0
          rms_threshold = 0.5
      [[slitedges]]
          edge_thresh = 50.0
          sync_predict = nearest
  [scienceframe]
      exprng = 90, None,
      [[process]]
          mask_cr = True
          sigclip = 5.0
          objlim = 2.0
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[findobj]]
          snr_thresh = 5.0
      [[skysub]]
          sky_sigrej = 5.0
          global_sky_std = False
          no_poly = True
  [sensfunc]
      polyorder = 7

.. _instr_par-gemini_flamingos1:

GEMINI-S FLAMINGOS (``gemini_flamingos1``)
------------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = gemini_flamingos1
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 20, None,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 1, 50,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          exprng = None, 60,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[wavelengths]]
          method = full_template
          lamps = ArI, ArII, ThAr, NeI,
          sigdetect = 3
          fwhm = 20
          reid_arxiv = magellan_fire_long.fits
          rms_threshold = 1.0
          match_toler = 5.0
      [[slitedges]]
          trace_thresh = 5.0
          sync_predict = nearest
      [[tilts]]
          tracethresh = 5
  [scienceframe]
      exprng = 20, None,
      [[process]]
          mask_cr = True
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[findobj]]
          snr_thresh = 5.0
          find_trim_edge = 50, 50,

.. _instr_par-gemini_flamingos2:

GEMINI-S FLAMINGOS (``gemini_flamingos2``)
------------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = gemini_flamingos2
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 20, None,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 50, None,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          exprng = 50, None,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          exprng = None, 30,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[wavelengths]]
          lamps = OH_NIRES,
          fwhm = 5
          rms_threshold = 0.5
          match_toler = 5.0
      [[slitedges]]
          edge_thresh = 200.0
          fit_min_spec_length = 0.4
          trace_thresh = 10.0
          sync_predict = nearest
      [[tilts]]
          tracethresh = 5
          spat_order = 4
  [scienceframe]
      exprng = 20, None,
      [[process]]
          mask_cr = True
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[findobj]]
          snr_thresh = 5.0
          find_trim_edge = 10, 10,
      [[skysub]]
          sky_sigrej = 5.0
  [sensfunc]
      algorithm = IR
      polyorder = 8
      [[IR]]
          telgridfile = TelFit_LasCampanas_3100_26100_R20000.fits

.. _instr_par-gemini_gmos_north_e2v:

GEMINI-N GMOS-N (``gemini_gmos_north_e2v``)
-------------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = gemini_gmos_north_e2v
      detnum = (1, 2, 3),
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = CuI, ArI, ArII,
          rms_threshold = 0.4
          nsnippet = 1
      [[slitedges]]
          edge_thresh = 100.0
          follow_span = 80
          fit_order = 3
      [[tilts]]
          tracethresh = 10.0
  [scienceframe]
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar

.. _instr_par-gemini_gmos_north_ham:

GEMINI-N GMOS-N (``gemini_gmos_north_ham``)
-------------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = gemini_gmos_north_ham
      detnum = (1, 2, 3),
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = CuI, ArI, ArII,
          rms_threshold = 0.4
          nsnippet = 1
      [[slitedges]]
          edge_thresh = 100.0
          follow_span = 80
          fit_order = 3
      [[tilts]]
          tracethresh = 10.0
  [scienceframe]
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar

.. _instr_par-gemini_gmos_north_ham_ns:

GEMINI-N GMOS-N (``gemini_gmos_north_ham_ns``)
----------------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = gemini_gmos_north_ham_ns
      detnum = (1, 2, 3),
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = CuI, ArI, ArII,
          rms_threshold = 0.4
          nsnippet = 1
      [[slitedges]]
          edge_thresh = 100.0
          follow_span = 80
          fit_order = 3
      [[tilts]]
          tracethresh = 10.0
  [scienceframe]
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar

.. _instr_par-gemini_gmos_south_ham:

GEMINI-S GMOS-S (``gemini_gmos_south_ham``)
-------------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = gemini_gmos_south_ham
      detnum = (1, 2, 3),
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = CuI, ArI, ArII,
          rms_threshold = 0.4
          nsnippet = 1
      [[slitedges]]
          edge_thresh = 100.0
          follow_span = 80
          fit_order = 3
          bound_detector = True
      [[tilts]]
          tracethresh = 10.0
  [scienceframe]
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar
  [sensfunc]
      algorithm = IR
      [[IR]]
          telgridfile = TelFit_LasCampanas_3100_26100_R20000.fits

.. _instr_par-gemini_gnirs_echelle:

GEMINI-N GNIRS (``gemini_gnirs_echelle``)
-----------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = gemini_gnirs_echelle
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = None, 30,
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = None, 30,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          exprng = None, 30,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[flatfield]]
          tweak_slits_thresh = 0.9
      [[tilts]]
          spat_order = 1
  [scienceframe]
      exprng = 30, None,
      [[process]]
          mask_cr = True
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[findobj]]
          find_trim_edge = 2, 2,
          maxnumber_sci = 2
          maxnumber_std = 1
      [[skysub]]
          bspline_spacing = 0.8
          global_sky_std = False
          no_poly = True
      [[extraction]]
          model_full_slit = True
  [sensfunc]
      algorithm = IR
      polyorder = 6
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-gemini_gnirs_ifu:

GEMINI-N GNIRS (``gemini_gnirs_ifu``)
-------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = gemini_gnirs_ifu
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = None, 30,
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = None, 30,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          exprng = None, 30,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[flatfield]]
          tweak_slits_thresh = 0.0
          tweak_slits_maxfrac = 0.0
          slit_trim = 2
          slit_illum_finecorr = False
      [[slitedges]]
          pad = 2
      [[tilts]]
          spat_order = 1
          spec_order = 1
  [scienceframe]
      exprng = 30, None,
      [[process]]
          mask_cr = True
          sigclip = 4.0
          objlim = 1.5
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[findobj]]
          find_trim_edge = 2, 2,
          maxnumber_sci = 2
          maxnumber_std = 1
      [[skysub]]
          global_sky_std = False
          no_poly = True
      [[extraction]]
          model_full_slit = True
          skip_extraction = True
      [[cube]]
          grating_corr = False
  [flexure]
      spec_maxshift = 0
  [sensfunc]
      algorithm = IR
      polyorder = 6
      [[UVIS]]
          extinct_correct = False
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-gtc_maat:

GTC OSIRIS (``gtc_maat``)
-------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = gtc_maat
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 180,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[flatfield]]
          slit_illum_finecorr = False
      [[wavelengths]]
          method = full_template
          lamps = XeI, HgI, NeI, ArI,
      [[slitedges]]
          sync_predict = nearest
          bound_detector = True
      [[tilts]]
          spat_order = 1
          spec_order = 1
  [scienceframe]
      exprng = 90, None,
      [[process]]
          mask_cr = True
          sigclip = 4.0
          objlim = 1.5
          use_biasimage = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[findobj]]
          maxnumber_std = 1
          skip_final_global = True
          skip_skysub = True
      [[skysub]]
          no_poly = True
      [[extraction]]
          skip_extraction = True
  [flexure]
      spec_maxshift = 3
  [sensfunc]
      [[UVIS]]
          extinct_correct = False

.. _instr_par-gtc_osiris:

GTC OSIRIS (``gtc_osiris``)
---------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = gtc_osiris
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 180,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = XeI, HgI, NeI, ArI,
      [[slitedges]]
          sync_predict = nearest
          bound_detector = True
      [[tilts]]
          spat_order = 5
          spec_order = 5
  [scienceframe]
      exprng = 90, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [reduce]
      [[findobj]]
          maxnumber_std = 1
      [[skysub]]
          no_poly = True

.. _instr_par-gtc_osiris_plus:

GTC OSIRIS (``gtc_osiris_plus``)
--------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = gtc_osiris_plus
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 180,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = XeI, HgI, NeI, ArI,
      [[slitedges]]
          sync_predict = nearest
          bound_detector = True
      [[tilts]]
          spat_order = 5
          spec_order = 5
  [scienceframe]
      exprng = 90, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [reduce]
      [[findobj]]
          maxnumber_std = 1
      [[skysub]]
          no_poly = True

.. _instr_par-jwst_nircam:

JWST NIRCAM (``jwst_nircam``)
-----------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = jwst_nircam
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          refframe = observed
  [scienceframe]
      [[process]]
          sigclip = 5.0
          objlim = 2.0
          noise_floor = 0.01
  [reduce]
      trim_edge = 0, 0,
      [[findobj]]
          find_trim_edge = 0, 0,
          maxnumber_sci = 2
          find_fwhm = 2.0
      [[skysub]]
          bspline_spacing = 1.2
          sky_sigrej = 4.0
          max_mask_frac = 0.95
      [[extraction]]
          boxcar_radius = 0.25
          sn_gauss = 6.0
          model_full_slit = True
          use_2dmodel_mask = False

.. _instr_par-jwst_nirspec:

JWST NIRSPEC (``jwst_nirspec``)
-------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = jwst_nirspec
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          refframe = observed
  [scienceframe]
      [[process]]
          sigclip = 5.0
          objlim = 2.0
          noise_floor = 0.01
  [reduce]
      trim_edge = 0, 0,
      [[findobj]]
          find_trim_edge = 0, 0,
          maxnumber_sci = 2
          find_fwhm = 2.0
      [[skysub]]
          bspline_spacing = 5.0
          sky_sigrej = 4.0
          mask_by_boxcar = True
          max_mask_frac = 0.95
      [[extraction]]
          boxcar_radius = 0.2
          sn_gauss = 5.0
          model_full_slit = True
          use_2dmodel_mask = False

.. _instr_par-keck_deimos:

KECK DEIMOS (``keck_deimos``)
-----------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_deimos
      detnum = (1, 5), (2, 6), (3, 7), (4, 8),
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              comb_sigrej = 10.0
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
      [[wavelengths]]
          lamps = ArI, NeI, KrI, XeI,
          match_toler = 2.5
          n_first = 3
      [[slitedges]]
          edge_thresh = 50.0
          follow_span = 1000
          fit_order = 3
          minimum_slit_length_sci = 4.0
          minimum_slit_gap = 0.25
      [[tilts]]
          tracethresh = 10
  [scienceframe]
      [[process]]
          mask_cr = True
          sigclip = 4.0
          objlim = 1.5
          use_biasimage = False
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar
  [sensfunc]
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-keck_esi:

KECK ESI (``keck_esi``)
-----------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_esi
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 1, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 300, None,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 60,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_sigrej = 3.0
          lamps = CuI, ArI, NeI, HgI, XeI, ArII,
          fwhm = 2.9
          fwhm_fromlines = True
          reid_arxiv = keck_esi_ECH.fits
          cc_thresh = 0.5
          cc_local_thresh = 0.5
          rms_threshold = 0.3
      [[slitedges]]
          edge_thresh = 5.0
          det_min_spec_length = 0.2
          max_shift_adj = 3.0
          fit_min_spec_length = 0.4
          left_right_pca = True
          pca_order = 3
          pca_sigrej = 1.5
          add_missed_orders = True
      [[tilts]]
          tracethresh = 10.0
  [scienceframe]
      exprng = 60, None,
      [[process]]
          satpix = nothing
          mask_cr = True
          sigclip = 20.0
          noise_floor = 0.01
  [reduce]
      [[findobj]]
          find_trim_edge = 4, 4,
          maxnumber_sci = 2
          maxnumber_std = 1
      [[extraction]]
          model_full_slit = True

.. _instr_par-keck_hires:

KECK HIRES (``keck_hires``)
---------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_hires
      detnum = (1, 2, 3),
  [calibrations]
      [[biasframe]]
          [[[process]]]
              overscan_method = median
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              overscan_method = median
              mask_cr = True
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              overscan_method = median
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              overscan_method = median
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              overscan_method = median
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              overscan_method = median
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              overscan_method = median
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
              use_pixelflat = False
              use_illumflat = False
      [[standardframe]]
          [[[process]]]
              overscan_method = median
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
              use_pixelflat = False
              use_illumflat = False
      [[flatfield]]
          tweak_slits_thresh = 0.9
      [[wavelengths]]
          method = echelle
          echelle = True
          ech_sigrej = 3.0
          lamps = ThAr,
          fwhm = 8.0
          cc_thresh = 0.5
          cc_local_thresh = 0.5
          rms_threshold = 0.5
          ech_separate_2d = True
      [[slitedges]]
          edge_thresh = 8.0
          max_shift_adj = 0.5
          fit_order = 8
          left_right_pca = True
          trace_thresh = 10.0
          max_nudge = 0.0
          dlength_range = 0.25
          length_range = 0.3
          overlap = True
      [[tilts]]
          tracethresh = 15
          spec_order = 5
  [scienceframe]
      [[process]]
          overscan_method = median
          mask_cr = True
          use_biasimage = False
          noise_floor = 0.01
          use_pixelflat = False
          use_illumflat = False
  [reduce]
      [[findobj]]
          find_trim_edge = 3, 3,
      [[skysub]]
          global_sky_std = False
      [[extraction]]
          model_full_slit = True
  [coadd1d]
      wave_method = log10
  [sensfunc]
      algorithm = IR
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-keck_kcrm:

KECK KCRM (``keck_kcrm``)
-------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_kcrm
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 0.01, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignment]]
          locations = 0.1, 0.3, 0.5, 0.7, 0.9,
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[flatfield]]
          spec_samp_coarse = 20.0
          tweak_slits_thresh = 0.0
          tweak_slits_maxfrac = 0.0
          slit_illum_relative = True
          slit_illum_ref_idx = 14
          slit_illum_smooth_npix = 5
          fit_2d_det_response = True
      [[wavelengths]]
          fwhm_spat_order = 2
      [[slitedges]]
          edge_thresh = 5
          fit_order = 4
          pad = 2
  [scienceframe]
      [[process]]
          mask_cr = True
          sigclip = 4.0
          objlim = 1.5
          use_biasimage = False
          noise_floor = 0.01
          use_specillum = True
  [reduce]
      [[skysub]]
          no_poly = True
      [[extraction]]
          skip_extraction = True
  [flexure]
      spec_maxshift = 3
  [sensfunc]
      [[UVIS]]
          extinct_correct = False

.. _instr_par-keck_kcwi:

KECK KCWI (``keck_kcwi``)
-------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_kcwi
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
              use_pattern = True
      [[darkframe]]
          exprng = 0.01, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
              use_pattern = True
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignment]]
          locations = 0.1, 0.3, 0.5, 0.7, 0.9,
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_illumflat = False
              use_pattern = True
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
              use_pattern = True
      [[flatfield]]
          spec_samp_coarse = 20.0
          tweak_slits_thresh = 0.0
          tweak_slits_maxfrac = 0.0
          slit_illum_relative = True
          slit_illum_ref_idx = 14
          slit_illum_smooth_npix = 5
          fit_2d_det_response = True
      [[wavelengths]]
          fwhm_spat_order = 2
      [[slitedges]]
          edge_thresh = 5
          fit_order = 4
          pad = 2
  [scienceframe]
      [[process]]
          mask_cr = True
          sigclip = 4.0
          objlim = 1.5
          use_biasimage = False
          noise_floor = 0.01
          use_specillum = True
          use_pattern = True
  [reduce]
      [[skysub]]
          no_poly = True
      [[extraction]]
          skip_extraction = True
  [flexure]
      spec_maxshift = 3
  [sensfunc]
      [[UVIS]]
          extinct_correct = False

.. _instr_par-keck_lris_blue:

KECK LRISb (``keck_lris_blue``)
-------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_lris_blue
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = None, 300,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = None, 300,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          exprng = None, 300,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = 1, 61,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
              spat_flexure_correct = True
      [[wavelengths]]
          method = full_template
          sigdetect = 10.0
          rms_threshold = 0.2
          n_first = 3
      [[slitedges]]
          edge_thresh = 15.0
          det_min_spec_length = 0.1
          fit_order = 3
          fit_min_spec_length = 0.2
          sync_center = gap
          minimum_slit_length = 3.0
          minimum_slit_length_sci = 5.0
  [scienceframe]
      exprng = 61, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01
          spat_flexure_correct = True
  [flexure]
      spec_method = boxcar
  [sensfunc]
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-keck_lris_blue_orig:

KECK LRISb (``keck_lris_blue_orig``)
------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_lris_blue_orig
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = None, 300,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = None, 300,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          exprng = None, 300,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = 1, 61,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
              spat_flexure_correct = True
      [[wavelengths]]
          method = full_template
          sigdetect = 10.0
          rms_threshold = 0.2
          n_first = 3
      [[slitedges]]
          edge_thresh = 15.0
          det_min_spec_length = 0.1
          fit_order = 3
          fit_min_spec_length = 0.2
          sync_center = gap
          minimum_slit_length = 3.0
          minimum_slit_length_sci = 5.0
  [scienceframe]
      exprng = 61, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01
          spat_flexure_correct = True
  [flexure]
      spec_method = boxcar
  [sensfunc]
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-keck_lris_red:

KECK LRISr (``keck_lris_red``)
------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_lris_red
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = None, 60,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = None, 60,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          exprng = None, 60,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = 1, 61,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
              spat_flexure_correct = True
      [[wavelengths]]
          sigdetect = 10.0
          rms_threshold = 0.2
      [[slitedges]]
          fit_order = 3
          sync_center = gap
          minimum_slit_length = 3.0
          minimum_slit_length_sci = 5.0
      [[tilts]]
          tracethresh = 25
          maxdev_tracefit = 1.0
          spat_order = 4
          spec_order = 7
          maxdev2d = 1.0
          sigrej2d = 5.0
  [scienceframe]
      exprng = 61, None,
      [[process]]
          mask_cr = True
          sigclip = 5.0
          objlim = 5.0
          noise_floor = 0.01
          spat_flexure_correct = True
  [reduce]
      [[skysub]]
          bspline_spacing = 0.8
  [flexure]
      spec_method = boxcar
  [sensfunc]
      algorithm = IR
      polyorder = 9
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-keck_lris_red_mark4:

KECK LRISr (``keck_lris_red_mark4``)
------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_lris_red_mark4
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = None, 60,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = None, 60,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          exprng = None, 60,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = 1, 61,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
              spat_flexure_correct = True
      [[wavelengths]]
          sigdetect = 10.0
          rms_threshold = 0.2
      [[slitedges]]
          fit_order = 3
          sync_center = gap
          minimum_slit_length = 3.0
          minimum_slit_length_sci = 5.0
      [[tilts]]
          tracethresh = 25
          maxdev_tracefit = 1.0
          spat_order = 4
          spec_order = 7
          maxdev2d = 1.0
          sigrej2d = 5.0
  [scienceframe]
      exprng = 61, None,
      [[process]]
          mask_cr = True
          sigclip = 5.0
          objlim = 5.0
          noise_floor = 0.01
          spat_flexure_correct = True
  [reduce]
      [[skysub]]
          bspline_spacing = 0.8
  [flexure]
      spec_method = boxcar
  [sensfunc]
      algorithm = IR
      polyorder = 9
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-keck_lris_red_orig:

KECK LRISr (``keck_lris_red_orig``)
-----------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_lris_red_orig
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = None, 60,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = None, 60,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          exprng = None, 60,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = 1, 61,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
              spat_flexure_correct = True
      [[wavelengths]]
          sigdetect = 10.0
          rms_threshold = 0.2
      [[slitedges]]
          fit_order = 3
          sync_center = gap
          minimum_slit_length = 3.0
          minimum_slit_length_sci = 5.0
      [[tilts]]
          tracethresh = 25
          maxdev_tracefit = 1.0
          spat_order = 4
          spec_order = 7
          maxdev2d = 1.0
          sigrej2d = 5.0
  [scienceframe]
      exprng = 61, None,
      [[process]]
          mask_cr = True
          sigclip = 5.0
          objlim = 5.0
          noise_floor = 0.01
          spat_flexure_correct = True
  [reduce]
      [[skysub]]
          bspline_spacing = 0.8
  [flexure]
      spec_method = boxcar
  [sensfunc]
      algorithm = IR
      polyorder = 9
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-keck_mosfire:

KECK MOSFIRE (``keck_mosfire``)
-------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_mosfire
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 1, None,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 1, None,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 20,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
      [[wavelengths]]
          lamps = OH_NIRES,
          fwhm = 5.0
          rms_threshold = 0.3
      [[slitedges]]
          edge_thresh = 50.0
          sync_predict = nearest
  [scienceframe]
      exprng = 20, None,
      [[process]]
          satpix = nothing
          mask_cr = True
          sigclip = 20.0
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
  [reduce]
      [[skysub]]
          bspline_spacing = 0.8
  [fluxcalib]
      extrap_sens = True
  [sensfunc]
      extrap_blu = 0.0
      extrap_red = 0.0
      algorithm = IR
      polyorder = 13
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-keck_nires:

KECK NIRES (``keck_nires``)
---------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_nires
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 61, None,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          exprng = 61, None,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          exprng = None, 60,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_norder_coeff = 6
          ech_sigrej = 3.0
          lamps = OH_NIRES,
          fwhm = 2.2
          fwhm_fromlines = True
          reid_arxiv = keck_nires.fits
          rms_threshold = 0.3
          n_final = 3, 4, 4, 4, 4,
      [[slitedges]]
          fit_min_spec_length = 0.4
          left_right_pca = True
          trace_thresh = 10.0
          fwhm_gaussian = 4.0
      [[tilts]]
          tracethresh = 10.0
  [scienceframe]
      exprng = 61, None,
      [[process]]
          satpix = nothing
          mask_cr = True
          sigclip = 20.0
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[skysub]]
          bspline_spacing = 0.8
      [[extraction]]
          boxcar_radius = 0.75
  [coadd1d]
      wave_method = log10
  [coadd2d]
      offsets = header
  [sensfunc]
      algorithm = IR
      polyorder = 8
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-keck_nirspec_low:

KECK NIRSPEC (``keck_nirspec_low``)
-----------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = keck_nirspec_low
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 20, None,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 20, None,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          exprng = None, 20,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[flatfield]]
          tweak_slits_thresh = 0.8
      [[wavelengths]]
          lamps = OH_NIRES,
          fwhm = 5.0
          rms_threshold = 0.2
      [[slitedges]]
          edge_thresh = 200.0
          sync_predict = nearest
  [scienceframe]
      exprng = 20, None,
      [[process]]
          satpix = nothing
          mask_cr = True
          sigclip = 20.0
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[skysub]]
          bspline_spacing = 0.8
  [sensfunc]
      algorithm = IR
      polyorder = 8
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-lbt_luci1:

LBT LUCI1 (``lbt_luci1``)
-------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = lbt_luci1
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[wavelengths]]
          lamps = OH_NIRES,
          fwhm = 5.0
          rms_threshold = 0.2
      [[slitedges]]
          edge_thresh = 300.0
          sync_predict = nearest
  [scienceframe]
      [[process]]
          satpix = nothing
          mask_cr = True
          sigclip = 20.0
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[skysub]]
          bspline_spacing = 0.8
      [[extraction]]
          std_prof_nsigma = 100.0

.. _instr_par-lbt_luci2:

LBT LUCI2 (``lbt_luci2``)
-------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = lbt_luci2
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[wavelengths]]
          lamps = OH_NIRES,
          fwhm = 5.0
          rms_threshold = 0.2
      [[slitedges]]
          edge_thresh = 300
          fit_order = 8
          sync_predict = nearest
  [scienceframe]
      [[process]]
          satpix = nothing
          mask_cr = True
          sigclip = 20.0
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[skysub]]
          bspline_spacing = 0.8
          global_sky_std = False
      [[extraction]]
          std_prof_nsigma = 100.0
          model_full_slit = True

.. _instr_par-lbt_mods1b:

LBT MODS1B (``lbt_mods1b``)
---------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = lbt_mods1b
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = 0, None,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = 0, None,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = 1, 200,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          lamps = XeI, KrI, ArI, HgI,
          sigdetect = 10.0
          rms_threshold = 0.4
      [[slitedges]]
          edge_thresh = 100.0
          sync_predict = nearest
      [[tilts]]
          maxdev_tracefit = 0.02
          spat_order = 5
          spec_order = 5
          maxdev2d = 0.02
  [scienceframe]
      exprng = 200, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar

.. _instr_par-lbt_mods1r:

LBT MODS1R (``lbt_mods1r``)
---------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = lbt_mods1r
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = 0, None,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = 0, None,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = 1, 200,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          lamps = ArI, NeI, KrI, XeI,
          fwhm = 10.0
          rms_threshold = 0.4
          match_toler = 2.5
          n_first = 3
      [[slitedges]]
          edge_thresh = 100.0
          sync_predict = nearest
      [[tilts]]
          maxdev_tracefit = 0.02
          spat_order = 5
          spec_order = 5
          maxdev2d = 0.02
  [scienceframe]
      exprng = 200, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar

.. _instr_par-lbt_mods2b:

LBT MODS2B (``lbt_mods2b``)
---------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = lbt_mods2b
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = 0, None,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = 0, None,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = 1, 200,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          lamps = XeI, KrI, ArI, HgI,
          sigdetect = 10.0
          rms_threshold = 0.4
      [[slitedges]]
          edge_thresh = 100.0
          sync_predict = nearest
      [[tilts]]
          maxdev_tracefit = 0.02
          spat_order = 5
          spec_order = 5
          maxdev2d = 0.02
  [scienceframe]
      exprng = 200, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar

.. _instr_par-lbt_mods2r:

LBT MODS2R (``lbt_mods2r``)
---------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = lbt_mods2r
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = 0, None,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = 0, None,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = 1, 200,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          lamps = ArI, NeI, KrI, XeI,
          fwhm = 10.0
          rms_threshold = 1.0
          match_toler = 2.5
          n_first = 3
      [[slitedges]]
          edge_thresh = 300.0
          sync_predict = nearest
      [[tilts]]
          maxdev_tracefit = 0.02
          spat_order = 5
          spec_order = 5
          maxdev2d = 0.02
  [scienceframe]
      exprng = 200, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar

.. _instr_par-ldt_deveny:

LDT DeVeny (``ldt_deveny``)
---------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = ldt_deveny
  [calibrations]
      bpm_usebias = True
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
              use_illumflat = False
      [[flatfield]]
          spec_samp_fine = 30
          pixelflat_min_wave = 3000.0
          slit_illum_finecorr = False
      [[wavelengths]]
          method = full_template
          lamps = use_header,
          fwhm_fromlines = True
          n_first = 3
          n_final = 5
          nsnippet = 1
      [[slitedges]]
          sync_predict = nearest
          bound_detector = True
          minimum_slit_length = 90.0
      [[tilts]]
          spat_order = 4
          spec_order = 5
  [scienceframe]
      [[process]]
          mask_cr = True
          sigclip = 5.0
          objlim = 2.0
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[findobj]]
          trace_npoly = 3
          snr_thresh = 50.0
          maxnumber_sci = 5
          maxnumber_std = 1
          find_fwhm = 3.5
      [[skysub]]
          sky_sigrej = 4.0
      [[extraction]]
          boxcar_radius = 1.8
          use_2dmodel_mask = False
  [flexure]
      spec_method = boxcar
      spec_maxshift = 30
  [sensfunc]
      [[UVIS]]
          polycorrect = False
          nresln = 15

.. _instr_par-magellan_fire:

MAGELLAN FIRE (``magellan_fire``)
---------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = magellan_fire
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 20, None,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 20, None,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          exprng = None, 60,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_norder_coeff = 6
          ech_sigrej = 3.0
          lamps = OH_FIRE_Echelle,
          sigdetect = 5, 10, 10, 10, 10, 20, 30, 30, 30, 30, 30, 10, 30, 30, 60, 30, 30, 10, 20, 30, 10,
          reid_arxiv = magellan_fire_echelle.fits
          cc_thresh = 0.35
          rms_threshold = 1.0
          match_toler = 30.0
          n_final = 3, 3, 3, 2, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 6, 6, 4,
      [[slitedges]]
          edge_thresh = 3.0
          max_shift_adj = 0.5
          fit_min_spec_length = 0.5
          left_right_pca = True
          pca_order = 3
          trace_thresh = 10.0
      [[tilts]]
          tracethresh = 5
  [scienceframe]
      exprng = 20, None,
      [[process]]
          satpix = nothing
          mask_cr = True
          sigclip = 20.0
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[findobj]]
          maxnumber_sci = 2
          maxnumber_std = 1
      [[extraction]]
          model_full_slit = True
  [coadd1d]
      wave_method = log10
  [sensfunc]
      algorithm = IR
      [[IR]]
          telgridfile = TelFit_LasCampanas_3100_26100_R20000.fits

.. _instr_par-magellan_fire_long:

MAGELLAN FIRE (``magellan_fire_long``)
--------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = magellan_fire_long
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 20, None,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 1, 50,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          exprng = None, 60,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[wavelengths]]
          method = full_template
          lamps = ArI, ArII, ThAr, NeI,
          sigdetect = 3
          fwhm = 20
          reid_arxiv = magellan_fire_long.fits
          rms_threshold = 1.0
          match_toler = 5.0
      [[slitedges]]
          trace_thresh = 10.0
          sync_predict = nearest
      [[tilts]]
          tracethresh = 5
  [scienceframe]
      exprng = 20, None,
      [[process]]
          mask_cr = True
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[findobj]]
          snr_thresh = 5
          find_trim_edge = 50, 50,
  [sensfunc]
      [[IR]]
          telgridfile = TelFit_LasCampanas_3100_26100_R20000.fits

.. _instr_par-magellan_mage:

MAGELLAN MagE (``magellan_mage``)
---------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = magellan_mage
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 20, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 20, None,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 20,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_sigrej = 3.0
          lamps = ThAr_MagE,
          fwhm = 3.0
          fwhm_fromlines = True
          reid_arxiv = magellan_mage.fits
          cc_thresh = 0.5
          cc_local_thresh = 0.5
          rms_threshold = 0.4
      [[slitedges]]
          edge_thresh = 10.0
          max_shift_adj = 3.0
          fit_min_spec_length = 0.3
          left_right_pca = True
      [[tilts]]
          tracethresh = 10.0
  [scienceframe]
      exprng = 20, None,
      [[process]]
          satpix = nothing
          mask_cr = True
          sigclip = 20.0
          noise_floor = 0.01
  [reduce]
      [[findobj]]
          find_trim_edge = 4, 4,
          maxnumber_sci = 2
          maxnumber_std = 1
      [[extraction]]
          model_full_slit = True
  [coadd1d]
      wave_method = log10

.. _instr_par-mdm_modspec:

HILTNER Echelle (``mdm_modspec``)
---------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = mdm_modspec
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              overscan_method = median
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              n_lohi = 1, 1,
              comb_sigrej = 3.0
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[flatfield]]
          slit_illum_finecorr = False
      [[wavelengths]]
          method = full_template
          lamps = ArI, XeI, NeI,
          reid_arxiv = mdm_modspec_1200_5100.fits
          n_final = 9
      [[slitedges]]
          sync_predict = nearest
          bound_detector = True
  [scienceframe]
      exprng = 10, 600,
      [[process]]
          mask_cr = True
          noise_floor = 0.01

.. _instr_par-mdm_osmos_mdm4k:

HILTNER MDM4K (``mdm_osmos_mdm4k``)
-----------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = mdm_osmos_mdm4k
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 120,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = ArI, XeI,
          sigdetect = 10.0
          reid_arxiv = mdm_osmos_mdm4k.fits
      [[slitedges]]
          sync_predict = nearest
  [scienceframe]
      exprng = 90, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01

.. _instr_par-mmt_binospec:

MMT BINOSPEC (``mmt_binospec``)
-------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = mmt_binospec
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 20, None,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 20, None,
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 100,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = HeI, NeI, ArI, ArII,
          fwhm = 5.0
          rms_threshold = 0.5
      [[slitedges]]
          sync_predict = nearest
      [[tilts]]
          tracethresh = 10.0
          spat_order = 6
          spec_order = 6
  [scienceframe]
      exprng = 20, None,
      [[process]]
          mask_cr = True
          sigclip = 5.0
          objlim = 2.0
          use_biasimage = False
          noise_floor = 0.01
  [reduce]
      [[skysub]]
          bspline_spacing = 0.8
          global_sky_std = False
  [flexure]
      spec_method = boxcar
  [sensfunc]
      polyorder = 7
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-mmt_bluechannel:

MMT Blue_Channel (``mmt_bluechannel``)
--------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = mmt_bluechannel
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 300, None,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 1, None,
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = None, 600,
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = None, 600,
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          exprng = 1, None,
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          exprng = None, 600,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
              use_illumflat = False
      [[wavelengths]]
          lamps = use_header,
          fwhm_fromlines = True
          rms_threshold = 0.5
      [[slitedges]]
          sync_predict = nearest
          bound_detector = True
  [scienceframe]
      [[process]]
          mask_cr = True
          sigclip = 5.0
          objlim = 2.0
          use_biasimage = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[skysub]]
          bspline_spacing = 0.8
          global_sky_std = False
  [sensfunc]
      polyorder = 7

.. _instr_par-mmt_mmirs:

MMT MMIRS (``mmt_mmirs``)
-------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = mmt_mmirs
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 30, None,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 60, None,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          exprng = 60, None,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          exprng = None, 60,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[wavelengths]]
          lamps = OH_NIRES,
          fwhm = 5
          rms_threshold = 0.5
          match_toler = 5.0
      [[slitedges]]
          edge_thresh = 100.0
          fit_min_spec_length = 0.4
          trace_thresh = 10.0
          sync_predict = nearest
          bound_detector = True
      [[tilts]]
          tracethresh = 5
          spat_order = 7
          spec_order = 5
  [scienceframe]
      exprng = 30, None,
      [[process]]
          mask_cr = True
          grow = 0.5
          sigclip = 5.0
          objlim = 2.0
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
  [reduce]
      [[findobj]]
          snr_thresh = 5.0
      [[skysub]]
          sky_sigrej = 5.0
  [sensfunc]
      algorithm = IR
      polyorder = 8
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-not_alfosc:

NOT ALFOSC (``not_alfosc``)
---------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = not_alfosc
  [calibrations]
      [[biasframe]]
          exprng = None, 1,
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
          [[[process]]]
              use_overscan = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_overscan = False
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 120,
          [[[process]]]
              mask_cr = True
              use_overscan = False
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = HeI, NeI, ArI,
          sigdetect = 10.0
      [[slitedges]]
          edge_thresh = 30
          sync_predict = nearest
          bound_detector = True
          minimum_slit_gap = 15
  [scienceframe]
      exprng = 10, None,
      [[process]]
          mask_cr = True
          use_overscan = False
          noise_floor = 0.01

.. _instr_par-not_alfosc_vert:

NOT ALFOSC (``not_alfosc_vert``)
--------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = not_alfosc_vert
  [calibrations]
      [[biasframe]]
          exprng = None, 1,
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
          [[[process]]]
              use_overscan = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_overscan = False
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 120,
          [[[process]]]
              mask_cr = True
              use_overscan = False
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = HeI, NeI, ArI,
          sigdetect = 10.0
      [[slitedges]]
          edge_thresh = 30
          sync_predict = nearest
          bound_detector = True
          minimum_slit_gap = 15
  [scienceframe]
      exprng = 10, None,
      [[process]]
          mask_cr = True
          use_overscan = False
          noise_floor = 0.01

.. _instr_par-ntt_efosc2:

NTT EFOSC2 (``ntt_efosc2``)
---------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = ntt_efosc2
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[flatfield]]
          tweak_slits_thresh = 0.9
      [[wavelengths]]
          method = full_template
          lamps = HeI, ArI,
          sigdetect = 10.0
          rms_threshold = 0.25
      [[slitedges]]
          edge_thresh = 75.0
          sync_predict = nearest
      [[tilts]]
          tracethresh = 25.0
  [scienceframe]
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [reduce]
      [[skysub]]
          sky_sigrej = 5.0
          global_sky_std = False
          no_poly = True
  [flexure]
      spec_method = boxcar

.. _instr_par-p200_dbsp_blue:

P200 DBSPb (``p200_dbsp_blue``)
-------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = p200_dbsp_blue
  [calibrations]
      bpm_usebias = True
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = None, 120,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 120,
          [[[process]]]
              combine = median
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = FeI, ArI, ArII,
      [[slitedges]]
          fit_min_spec_length = 0.55
          sync_predict = nearest
  [scienceframe]
      exprng = 90, None,
      [[process]]
          combine = median
          mask_cr = True
          noise_floor = 0.01
  [sensfunc]
      [[UVIS]]
          nresln = 5

.. _instr_par-p200_dbsp_red:

P200 DBSPr (``p200_dbsp_red``)
------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = p200_dbsp_red
  [calibrations]
      bpm_usebias = True
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = None, 120,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 120,
          [[[process]]]
              combine = median
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = ArI, ArII, NeI, HeI,
      [[slitedges]]
          sync_predict = nearest
  [scienceframe]
      exprng = 90, None,
      [[process]]
          combine = median
          mask_cr = True
          sigclip = 4.0
          objlim = 1.5
          noise_floor = 0.01
  [sensfunc]
      [[UVIS]]
          polycorrect = False
      [[IR]]
          telgridfile = TelFit_Lick_3100_11100_R10000.fits

.. _instr_par-p200_tspec:

P200 TSPEC (``p200_tspec``)
---------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = p200_tspec
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 0, None,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 100, None,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          exprng = 100, None,
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          exprng = None, 60,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_norder_coeff = 6
          ech_sigrej = 3.0
          lamps = OH_NIRES,
          fwhm = 2.9
          fwhm_fromlines = True
          reid_arxiv = p200_triplespec.fits
          rms_threshold = 0.3
          n_final = 3, 4, 4, 4, 4,
      [[slitedges]]
          fit_min_spec_length = 0.3
          left_right_pca = True
          trace_thresh = 5.0
          fwhm_gaussian = 4.0
      [[tilts]]
          tracethresh = 10.0
  [scienceframe]
      exprng = 60, None,
      [[process]]
          satpix = nothing
          mask_cr = True
          sigclip = 20.0
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[findobj]]
          maxnumber_sci = 2
          maxnumber_std = 1
      [[skysub]]
          bspline_spacing = 0.8
      [[extraction]]
          boxcar_radius = 0.75
          model_full_slit = True
  [coadd1d]
      wave_method = log10
  [sensfunc]
      algorithm = IR
      polyorder = 8
      [[IR]]
          telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

.. _instr_par-shane_kast_blue:

SHANE KASTb (``shane_kast_blue``)
---------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = shane_kast_blue
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = None, 61,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = 0, None,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = 0, None,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = 1, 61,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = CdI, HgI, HeI,
          rms_threshold = 0.2
          match_toler = 2.5
          n_first = 3
      [[slitedges]]
          sync_predict = nearest
          bound_detector = True
      [[tilts]]
          maxdev_tracefit = 0.02
          spec_order = 5
          maxdev2d = 0.02
  [scienceframe]
      exprng = 61, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar
      spectrum = sky_kastb_600.fits
  [sensfunc]
      [[IR]]
          telgridfile = TelFit_Lick_3100_11100_R10000.fits

.. _instr_par-shane_kast_red:

SHANE KASTr (``shane_kast_red``)
--------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = shane_kast_red
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = None, 61,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = 0, None,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = 0, None,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = 1, 61,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          lamps = NeI, HgI, HeI, ArI,
      [[slitedges]]
          sync_predict = nearest
          bound_detector = True
  [scienceframe]
      exprng = 61, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar
  [sensfunc]
      [[IR]]
          telgridfile = TelFit_Lick_3100_11100_R10000.fits

.. _instr_par-shane_kast_red_ret:

SHANE KASTr (``shane_kast_red_ret``)
------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = shane_kast_red_ret
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = None, 61,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          exprng = 0, None,
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          exprng = 0, None,
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          exprng = 1, 61,
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[wavelengths]]
          lamps = NeI, HgI, HeI, ArI,
          rms_threshold = 0.2
          use_instr_flag = True
      [[slitedges]]
          sync_predict = nearest
          bound_detector = True
  [scienceframe]
      exprng = 61, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar
  [sensfunc]
      [[IR]]
          telgridfile = TelFit_Lick_3100_11100_R10000.fits

.. _instr_par-soar_goodman_blue:

SOAR blue (``soar_goodman_blue``)
---------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = soar_goodman_blue
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = None, 30,
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 120,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
      [[wavelengths]]
          lamps = NeI, ArI, HgI,
          fwhm = 5.0
          rms_threshold = 0.5
      [[slitedges]]
          sync_predict = nearest
          bound_detector = True
  [scienceframe]
      exprng = 90, None,
      [[process]]
          mask_cr = True
          use_biasimage = False
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar
  [sensfunc]
      [[IR]]
          telgridfile = TelFit_LasCampanas_3100_26100_R20000.fits

.. _instr_par-soar_goodman_red:

SOAR red (``soar_goodman_red``)
-------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = soar_goodman_red
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = None, 30,
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 120,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
      [[flatfield]]
          slit_illum_finecorr = False
      [[wavelengths]]
          lamps = NeI, ArI, HgI,
          fwhm = 5.0
          rms_threshold = 0.5
      [[slitedges]]
          sync_predict = nearest
          bound_detector = True
  [scienceframe]
      exprng = 90, None,
      [[process]]
          mask_cr = True
          use_biasimage = False
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar
  [sensfunc]
      [[IR]]
          telgridfile = TelFit_LasCampanas_3100_26100_R20000.fits

.. _instr_par-tng_dolores:

TNG DOLORES (``tng_dolores``)
-----------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = tng_dolores
  [calibrations]
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[tiltframe]]
          [[[process]]]
              clip = False
              use_pixelflat = False
              use_illumflat = False
              subtract_continuum = True
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              mask_cr = True
              noise_floor = 0.01
      [[slitedges]]
          sync_predict = nearest
  [scienceframe]
      exprng = 1, None,
      [[process]]
          mask_cr = True
          noise_floor = 0.01

.. _instr_par-vlt_fors2:

VLT FORS2 (``vlt_fors2``)
-------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = vlt_fors2
  [calibrations]
      [[biasframe]]
          [[[process]]]
              overscan_method = median
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              overscan_method = median
              mask_cr = True
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              overscan_method = median
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              overscan_method = median
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              overscan_method = median
      [[alignframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              overscan_method = median
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              overscan_method = median
              mask_cr = True
              noise_floor = 0.01
      [[standardframe]]
          [[[process]]]
              overscan_method = median
              mask_cr = True
              noise_floor = 0.01
      [[flatfield]]
          tweak_slits_thresh = 0.9
      [[wavelengths]]
          lamps = HeI, ArI,
          sigdetect = 10.0
          rms_threshold = 0.25
      [[slitedges]]
          edge_thresh = 50.0
          max_shift_adj = 0.5
          fit_order = 3
      [[tilts]]
          tracethresh = 25.0
  [scienceframe]
      [[process]]
          mask_cr = True
          noise_floor = 0.01
  [flexure]
      spec_method = boxcar
  [sensfunc]
      algorithm = IR
      [[IR]]
          telgridfile = TelFit_Paranal_VIS_9800_25000_R25000.fits

.. _instr_par-vlt_sinfoni:

VLT SINFONI (``vlt_sinfoni``)
-----------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = vlt_sinfoni
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 20, None,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = 20, None,
          [[[process]]]
              mask_cr = True
              sigclip = 20.0
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              mask_cr = True
              sigclip = 20.0
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              sigclip = 20.0
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 20,
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = OH_FIRE_Echelle,
          fwhm = 5.0
          reid_arxiv = vlt_sinfoni_K.fits
          rms_threshold = 0.3
          nsnippet = 1
      [[slitedges]]
          edge_thresh = 50.0
          sync_predict = nearest
          rm_slits = 1:1024:983
      [[tilts]]
          tracethresh = 5.0
  [scienceframe]
      exprng = 20, None,
      [[process]]
          satpix = nothing
          mask_cr = True
          sigclip = 20.0
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
  [reduce]
      [[findobj]]
          find_fwhm = 10
          skip_second_find = True
      [[skysub]]
          bspline_spacing = 0.9
          global_sky_std = False
      [[extraction]]
          sn_gauss = 5.0
          model_full_slit = True
  [sensfunc]
      algorithm = IR
      polyorder = 7
      [[IR]]
          telgridfile = TelFit_Paranal_NIR_9800_25000_R25000.fits

.. _instr_par-vlt_xshooter_nir:

VLT XShooter_NIR (``vlt_xshooter_nir``)
---------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = vlt_xshooter_nir
  [calibrations]
      [[biasframe]]
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_biasimage = False
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[standardframe]]
          [[[process]]]
              use_biasimage = False
              use_overscan = False
              noise_floor = 0.01
              use_illumflat = False
      [[flatfield]]
          tweak_slits_thresh = 0.9
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_nspec_coeff = 5
          ech_norder_coeff = 5
          ech_sigrej = 3.0
          lamps = OH_XSHOOTER,
          sigdetect = 10.0
          fwhm = 5.0
          reid_arxiv = vlt_xshooter_nir.fits
          cc_thresh = 0.5
          cc_local_thresh = 0.5
          rms_threshold = 0.6
          qa_log = False
      [[slitedges]]
          edge_thresh = 50.0
          max_shift_adj = 0.5
          fit_order = 8
          fit_min_spec_length = 0.5
          left_right_pca = True
          trace_thresh = 10.0
          length_range = 0.3
      [[tilts]]
          tracethresh = 25.0
          maxdev_tracefit = 0.04
          maxdev2d = 0.04
          rm_continuum = True
  [scienceframe]
      [[process]]
          satpix = nothing
          mask_cr = True
          sigclip = 20.0
          use_biasimage = False
          use_overscan = False
          noise_floor = 0.01
          use_illumflat = False
  [reduce]
      [[findobj]]
          trace_npoly = 8
          maxnumber_sci = 2
          maxnumber_std = 1
      [[skysub]]
          bspline_spacing = 0.8
          global_sky_std = False
      [[extraction]]
          model_full_slit = True
  [coadd1d]
      wave_method = log10
  [sensfunc]
      algorithm = IR
      polyorder = 8
      [[IR]]
          telgridfile = TelFit_Paranal_NIR_9800_25000_R25000.fits

.. _instr_par-vlt_xshooter_uvb:

VLT XShooter_UVB (``vlt_xshooter_uvb``)
---------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = vlt_xshooter_uvb
  [calibrations]
      [[biasframe]]
          [[[process]]]
              overscan_method = median
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              overscan_method = median
              mask_cr = True
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              overscan_method = median
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              overscan_method = median
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              overscan_method = median
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              overscan_method = median
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              overscan_method = median
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
              use_pixelflat = False
              use_illumflat = False
      [[standardframe]]
          [[[process]]]
              overscan_method = median
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
      [[flatfield]]
          tweak_slits_thresh = 0.9
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_sigrej = 3.0
          lamps = ThAr_XSHOOTER_UVB,
          sigdetect = 3.0
          fwhm = 3.8
          fwhm_fromlines = True
          reid_arxiv = vlt_xshooter_uvb1x1.fits
          cc_thresh = 0.5
          cc_local_thresh = 0.5
          rms_threshold = 0.7
          n_final = 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      [[slitedges]]
          edge_thresh = 8.0
          max_shift_adj = 0.5
          left_right_pca = True
          trace_thresh = 10.0
          length_range = 0.3
  [scienceframe]
      [[process]]
          overscan_method = median
          mask_cr = True
          noise_floor = 0.01
  [reduce]
      [[findobj]]
          find_trim_edge = 3, 3,
          maxnumber_sci = 2
          maxnumber_std = 1
      [[skysub]]
          bspline_spacing = 0.5
          global_sky_std = False
      [[extraction]]
          model_full_slit = True
  [coadd1d]
      wave_method = log10
  [sensfunc]
      algorithm = IR
      polyorder = 8
      [[IR]]
          telgridfile = TelFit_LasCampanas_3100_26100_R20000.fits

.. _instr_par-vlt_xshooter_vis:

VLT XShooter_VIS (``vlt_xshooter_vis``)
---------------------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = vlt_xshooter_vis
  [calibrations]
      [[biasframe]]
          [[[process]]]
              overscan_method = median
              combine = median
              use_biasimage = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          [[[process]]]
              overscan_method = median
              mask_cr = True
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          [[[process]]]
              overscan_method = median
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              overscan_method = median
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          [[[process]]]
              overscan_method = median
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[alignframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              overscan_method = median
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              overscan_method = median
              satpix = nothing
              use_biasimage = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              overscan_method = median
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
              use_pixelflat = False
              use_illumflat = False
      [[standardframe]]
          [[[process]]]
              overscan_method = median
              mask_cr = True
              use_biasimage = False
              noise_floor = 0.01
      [[flatfield]]
          tweak_slits_thresh = 0.9
      [[wavelengths]]
          method = reidentify
          echelle = True
          ech_sigrej = 3.0
          lamps = ThAr_XSHOOTER_VIS,
          fwhm = 8.0
          fwhm_fromlines = True
          reid_arxiv = vlt_xshooter_vis1x1.fits
          cc_thresh = 0.5
          cc_local_thresh = 0.5
          rms_threshold = 1.2
          n_final = 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3,
      [[slitedges]]
          edge_thresh = 8.0
          max_shift_adj = 0.5
          fit_order = 8
          left_right_pca = True
          trace_thresh = 10.0
          length_range = 0.3
      [[tilts]]
          tracethresh = 15
          spec_order = 5
  [scienceframe]
      [[process]]
          overscan_method = median
          mask_cr = True
          noise_floor = 0.01
  [reduce]
      [[findobj]]
          find_trim_edge = 3, 3,
          maxnumber_sci = 2
          maxnumber_std = 1
      [[skysub]]
          bspline_spacing = 0.5
          global_sky_std = False
      [[extraction]]
          model_full_slit = True
  [coadd1d]
      wave_method = log10
  [sensfunc]
      algorithm = IR
      polyorder = 8
      [[IR]]
          telgridfile = TelFit_Paranal_VIS_4900_11100_R25000.fits

.. _instr_par-wht_isis_blue:

WHT ISISb (``wht_isis_blue``)
-----------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = wht_isis_blue
  [calibrations]
      bpm_usebias = True
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = None, 120,
          [[[process]]]
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
          [[[process]]]
              use_overscan = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_overscan = False
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 120,
          [[[process]]]
              mask_cr = True
              use_overscan = False
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = NeI, ArI, ArII, CuI,
          sigdetect = 10.0
          n_first = 3
          n_final = 5
      [[slitedges]]
          sync_predict = nearest
  [scienceframe]
      exprng = 90, None,
      [[process]]
          mask_cr = True
          use_overscan = False
          noise_floor = 0.01

.. _instr_par-wht_isis_red:

WHT ISISr (``wht_isis_red``)
----------------------------
Alterations to the default parameters are:

.. code-block:: ini

  [rdx]
      spectrograph = wht_isis_red
  [calibrations]
      bpm_usebias = True
      [[biasframe]]
          exprng = None, 0.001,
          [[[process]]]
              combine = median
              use_biasimage = False
              use_overscan = False
              shot_noise = False
              use_pixelflat = False
              use_illumflat = False
      [[darkframe]]
          exprng = 999999, None,
          [[[process]]]
              mask_cr = True
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[arcframe]]
          exprng = None, 120,
          [[[process]]]
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[tiltframe]]
          [[[process]]]
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pixelflatframe]]
          [[[process]]]
              combine = median
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[pinholeframe]]
          exprng = 999999, None,
          [[[process]]]
              use_overscan = False
      [[alignframe]]
          [[[process]]]
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[traceframe]]
          [[[process]]]
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[illumflatframe]]
          [[[process]]]
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[lampoffflatsframe]]
          [[[process]]]
              satpix = nothing
              use_overscan = False
              use_pixelflat = False
              use_illumflat = False
      [[skyframe]]
          [[[process]]]
              mask_cr = True
              use_overscan = False
              noise_floor = 0.01
      [[standardframe]]
          exprng = None, 120,
          [[[process]]]
              mask_cr = True
              use_overscan = False
              noise_floor = 0.01
      [[wavelengths]]
          method = full_template
          lamps = NeI, ArI, ArII, CuI,
          sigdetect = 10.0
      [[slitedges]]
          sync_predict = nearest
  [scienceframe]
      exprng = 90, None,
      [[process]]
          mask_cr = True
          use_overscan = False
          noise_floor = 0.01

