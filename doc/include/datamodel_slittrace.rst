Version: 1.1.1

=====================  =================  ==========  ===========================================================================================================================
Obj Key                Obj Type           Array Type  Description                                                                                                                                                                         
=====================  =================  ==========  ===========================================================================================================================
``spat_id``            ndarray            int         Slit ID number from SPAT measured at half way point.
``maskdef_id``         ndarray            int         Slit ID number from slitmask design (implemented only for :doc:`deimos`).
``left_init``          ndarray            float       Spatial coordinates (pixel indices) of all left edges, one per slit.
``right_init``         ndarray            float       Spatial coordinates (pixel indices) of all right edges, one per slit.
``left_tweak``         ndarray            float       Spatial coordinates (pixel indices) of all left edges, one per slit.  These traces have been adjusted by the flat-field.
``right_tweak``        ndarray            float       Spatial coordinates (pixel indices) of all right edges, one per slit.  These traces have been adjusted by the flat-field.
``center``             ndarray            float       Spatial coordinates of the slit centers from ``left_init`` and ``right_init``.
``mask_init``          ndarray            int         Bit mask for slits at instantiation. Used to reset.
``mask``               ndarray            int         Bit mask for slits (fully good slits have 0 value).
``specmin``            ndarray            float       Minimum spectral position allowed for each slit/order.
``specmax``            ndarray            float       Maximum spectral position allowed for each slit/order.
=====================  =================  ==========  ===========================================================================================================================

Additional `astropy.io.fits.BinTableHDU`_ for :doc:`deimos` reduction.

=====================  =================  ==========  =================================================================================
Obj Key                Obj Type           Array Type  Description
=====================  =================  ==========  =================================================================================
``TRACEID``            ndarray            int         Trace ID Number.
``TRACESROW``          ndarray            int         Spectral row for provided left and right edges.
``TRACELPIX``          ndarray            float       Spatial pixel coordinate for left edge.
``TRACERPIX``          ndarray            float       Spatial pixel coordinate for right edge.
``SLITID``             ndarray            int         Slit ID Number (``maskdef_id``).
``SLITLOPT``           ndarray            float       Left edge of the slit in pixel from optical model before x-correlation.
``SLITROPT``           ndarray            float       Right edge of the slit in pixel from optical model before x-correlation.
``SLITRA``             ndarray            float       Right ascension of the slit center (deg).
``SLITDEC``            ndarray            float       Declination of the slit center (deg).
``SLITLEN``            ndarray            float       Slit length (arcsec).
``SLITWID``            ndarray            float       Slit width (arcsec).
``SLITPA``             ndarray            float       Slit position angle on sky (deg from N through E).
``ALIGN``              ndarray            int         Slit used for alignment (1-yes; 0-no).
``OBJID``              ndarray            int         Object ID Number.
``OBJRA``              ndarray            float       Right ascension of the object (deg).
``OBJDEC``             ndarray            float       Declination of the object (deg).
=====================  =================  ==========  =================================================================================