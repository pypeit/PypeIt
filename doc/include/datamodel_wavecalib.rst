
Version 1.1.2

========================  ==============================  =========  ===================================================================================================================================
HDU Name                  HDU Type                        Data Type  Description                                                                                                                        
========================  ==============================  =========  ===================================================================================================================================
``PRIMARY``               `astropy.io.fits.PrimaryHDU`_   ...        Empty data HDU.  Contains basic header information.                                                                                
``SPAT_IDS``              `astropy.io.fits.ImageHDU`_     integer    Slit spat_ids. Named distinctly from that in WaveFit                                                                               
``SPAT_ID-?__WAVEFIT``    `astropy.io.fits.BinTableHDU`_  ...        :class:`~pypeit.core.wavecal.wv_fitting.WaveFit` result for ``slit_id=?``                                                          
``SPAT_ID-?__PYPEITFIT``  `astropy.io.fits.BinTableHDU`_  ...        :class:`~pypeit.core.fitting.PypeItFit` element of the :class:`~pypeit.core.wavecal.wv_fitting.WaveFit` datamodel for ``slit_id=?``
...                       ...                             ...        ...                                                                                                                                
``ARC_SPECTRA``           `astropy.io.fits.ImageHDU`_     floating   2D array: 1D extracted spectra, slit by slit (nspec, nslits)                                                                       
========================  ==============================  =========  ===================================================================================================================================
