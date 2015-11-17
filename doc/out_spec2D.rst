.. highlight:: rest

*************
Spec2D Output 
*************

.. index:: spec2d

During the data reduction proceess, PYPIT creates a series
of 2D spectral images prior to extraction of 1D spectra.
And, of course, several of these 2D images may have greater
value for analysis than the 1D spectra.  For each on-source
exposure, PYPIT outputs a series of these images, with the
number set by the :ref:`outputs-reduction-mode`.  
The following table describes the possible products:



============  ====================================  ===============
2D Spec Type  Description                           Written?
============  ====================================  ===============
ivar          Inverse variance image; sky+detector  standard, full
mask          Mask image                            full
objmodel      Model of the object(s) flux           full
processed     Bias-subtracted, flat-fielded image   full
residual      Residual image; data-model            full
sensitivity   Sensitivity image for fluxing         full
              (surface brightness)
skysub        Sky-subtracted, processed image       standard, full
skymodel      Model of sky emission on detector     full
waveimg       Wavelength image.  Vacuum, flexure,   standard, full
              and helio-centric corrected
============  ====================================  ===============

