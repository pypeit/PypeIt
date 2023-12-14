.. include:: ../include/links.rst

.. _slits:

=====
Slits
=====


Overview
========

This file describes the ``Slits`` object.
It contains the main information on the traced slit edges, organized into left-right slit pairs.


This is written to disk as a multi-extension FITS file prefixed by
``Slits`` in the ``Calibrations/`` folder.
See :ref:`calib-naming` for the naming convention.


Viewing
=======

The preferred way to view the slit edges information contained in ``Slits`` is as follows:

.. code-block:: python

    from astropy.io import fits
    from astropy.table import Table

    hdu = fits.open('Slits_A_1_DET01.fits')
    Table(hdu['SLITS'].data)

This will show a table that looks like this:

.. code-block:: console

    spat_id maskdef_id             left_init [4096]                        right_init [4096]                        left_tweak [4096]                        right_tweak [4096]                         center [4096]               mask_init  mask specmin specmax
    int64    int64                    float64                                  float64                                  float64                                  float64                                  float64                    int16   int16 float64 float64
    ------- ---------- ---------------------------------------- ---------------------------------------- ---------------------------------------- ---------------------------------------- ---------------------------------------- --------- ----- ------- -------
         31     699084                 5.0 .. 37.18657261586156 29.277545928955078 .. 62.018354415893555                 5.0 .. 37.18657261586156 29.277545928955078 .. 62.018354415893555   17.13877296447754 .. 49.60246351587756         0     0    -inf     inf
        111     699078   33.88330268859863 .. 65.59163284301758  135.11531257629395 .. 167.1435375213623     37.948677458954 .. 69.65700761337295  135.11531257629395 .. 167.1435375213623  84.49930763244629 .. 116.36758518218994         0     0    -inf     inf
        217     699086 139.54998016357422 .. 170.55170822143555  243.45840644836426 .. 274.8020210266113  143.1504450626825 .. 174.15217312054384  243.45840644836426 .. 274.8020210266113 191.50419330596924 .. 222.67686462402344         0     0    -inf     inf
        322     699091   248.21048545837402 .. 278.352352142334   345.8145122528076 .. 376.4430561065674   251.7797389884793 .. 281.9216056724392   345.8145122528076 .. 376.4430561065674   297.0124988555908 .. 327.3977041244507         0     0    -inf     inf
        457     699080  350.76593017578125 .. 379.9788398742676   514.4993572235107 .. 544.0507583618164  355.0206685112119 .. 384.23357820969824    513.3126446738094 .. 542.864045812115     432.632643699646 .. 462.014799118042         0     0    -inf     inf

In addition, if reducing data from these :ref:`slitmask_info_instruments`
and slit-mask design matching is performed (see e.g., :ref:`deimos-mask-matching`
for DEIMOS and :ref:`mosfire-edge-tracing` for MOSFIRE), a second
`astropy.io.fits.BinTableHDU`_ is written to disk.

.. code-block:: console

    Table(hdu['MASKDEF_DESIGNTAB'].data)

    TRACEID TRACESROW     TRACELPIX          TRACERPIX      SPAT_ID  MASKDEF_ID       SLITLOPT           SLITROPT         SLITRA      SLITDEC        SLITLEN       SLITWID  SLITPA ALIGN OBJID     OBJRA        OBJDEC   OBJNAME  OBJMAG OBJMAG_BAND OBJ_TOPDIST OBJ_BOTDIST
     int64    int64        float64            float64        int64     int64        float64            float64         float64      float64        float64       float64 float64 int16 int64    float64      float64    str32  float64    str32      float64     float64
    ------- --------- ------------------ ------------------ ------- ----------- ------------------- ------------------ ------------ ----------- ------------------ ------- ------- ----- ------ ------------ ----------- ------- ------- ----------- ----------- -----------
          0      2048                5.0  56.08821487426758      31      699084 -21.955958625056724  69.84514502686432 358.68126787  42.3347698 10.682000160217285     1.0     0.0     0 733790  358.6812625   42.334675 3003915   21.39           I         5.0       5.682
          1      2048 59.798553466796875 161.30295944213867     111      699078   72.84633999813946  175.5490274983032 358.69310956 42.33811871 11.979000091552734     1.0     0.0     0 733784 358.69310417 42.33803333 3003737   22.54           I       5.682       6.297
          2      2048  164.9419937133789 269.14466094970703     217      699086  178.11184062306785 283.34493633572674 358.67168891 42.34139707 12.293000221252441     1.0     0.0     0 733792 358.67168333 42.34143889 3104178   21.19           I       6.297       5.996
          3      2048  272.8930130004883  370.9837417602539     322      699091   286.1885530241663   384.710234625776 358.65265564 42.34461342 11.529000282287598     1.0     0.0     0 733797    358.65265 42.34467778 3104468   20.24           I       5.996       5.533
          4      2048 374.71178817749023    538.83056640625     457      699080  387.93013336193303  552.8797949320767 358.68945391 42.34933141 19.270000457763672     1.0     0.0     0 733786    358.68945 42.34819167 3103868   22.44           I       5.533      13.737

See :ref:`slittrace-datamodel` for a description of the columns.

.. _slittrace-datamodel:

Current SlitTrace Data Model
============================

Internally, the ``Slits`` object is held in :class:`~pypeit.slittrace.SlitTraceSet`,
which subclasses from :class:`~pypeit.datamodel.DataContainer`.

Here is its datamodel, which is written as an `astropy.io.fits.BinTableHDU`_.
In addition, if reducing :doc:`../spectrographs/deimos` or
:doc:`../spectrographs/mosfire` data and slit-mask design matching is performed,
another  `astropy.io.fits.BinTableHDU`_ is written to a fits extension named
`MASKDEF_DESIGNTAB`.

.. include:: ../include/datamodel_slittraceset.rst

