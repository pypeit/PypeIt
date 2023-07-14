
:orphan:

.. include:: include/links.rst

.. _history:

===============
History Keyword
===============

PypeIt will update the ``HISTORY`` keyword in the FITS header with
information about how a file was created and what raw data the file was
reduced from.

Each history entry will start with a date in ISO 8601 format and the 
string PypeIt. Because of the limited length allowed in the FITS
standard, entries may wrap to the next ``HISTORY`` keyword. To make it 
easier to detect if this has happened to a filename, all filenames are 
wrapped in double quotes ``"``.

Examples
--------

Reduction:

.. code-block:: console

  HISTORY 2021-03-05T23:56 PypeIt Reducing target HIP15339                        
  HISTORY Combining frames:                                                       
  HISTORY "S20161108S0087.fits.gz"                                                
  HISTORY "S20161108S0090.fits.gz"                                                
  HISTORY Subtracted background from frames:                                      
  HISTORY "S20161108S0088.fits.gz"                                                
  HISTORY "S20161108S0089.fits.gz"                                                
  HISTORY Callibration frames:                                                    
  HISTORY arc,science,tilt "S20161108S0069.fits.gz"                               
  HISTORY arc,science,tilt "S20161108S0070.fits.gz"                               
  HISTORY arc,science,tilt "S20161108S0071.fits.gz"                               
  HISTORY arc,science,tilt "S20161108S0072.fits.gz"                               
  HISTORY pixelflat,trace "S20161108S0078.fits.gz"                                

Coadding:

.. code-block:: console

  HISTORY 2021-01-23T02:12 PypeIt Coadded 4 objects from 3 spec1d files           
  HISTORY File 0 "spec1d_DE.20170425.53065-dra11_DEIMOS_2017Apr25T144418.240.fits"  
  HISTORY File 1 "spec1d_DE.20170425.51771-dra11_DEIMOS_2017Apr25T142245.350.fits"  
  HISTORY File 2 "spec1d_DE.20170425.50487-dra11_DEIMOS_2017Apr25T140121.014.fits"  
  HISTORY Object ID SPAT0692-SLIT0704-DET08 from file 0                           
  HISTORY Object ID SPAT0695-SLIT0706-DET04 from file 2                           
  HISTORY Object ID SPAT0691-SLIT0704-DET08 from file 2                           
  HISTORY Object ID SPAT0695-SLIT0706-DET04 from file 1

Fluxing:

.. code-block:: console

  HISTORY 2021-03-09T01:21 PypeIt Flux calibration "sens_b24-Feige66_KASTb_2015May
  HISTORY 20T041246.960.fits"                                                     

