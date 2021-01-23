.. include:: include/links.rst

.. _master_slits:

===========
MasterSlits
===========


Overview
========

This file describes the MasterSlits object.
It contains the main information on the traced slit edges, organized into left-right slit pairs.


See below for the `Current SlitTrace Data Model`_.
This is written to disk as a multi-extension FITS file prefixed by
MasterSlits in the Masters/ folder.
See the :ref:`master-naming` docs for more.


Viewing
=======

This is the preferred way to view the slit edges information contained in MasterSlits:

.. code-block:: python

    from astropy.io import fits
    from astropy.table import Table

    hdu=fits.open('MasterSlits_A_1_01.fits')
    Table(hdu[1].data)

This will show a table that looks like this:

.. code-block:: python

    spat_id maskdef_id             left_init [4096]                        right_init [4096]                        left_tweak [4096]                        right_tweak [4096]                         center [4096]               mask_init  mask specmin specmax
    int64    int64                    float64                                  float64                                  float64                                  float64                                  float64                    int16   int16 float64 float64
    ------- ---------- ---------------------------------------- ---------------------------------------- ---------------------------------------- ---------------------------------------- ---------------------------------------- --------- ----- ------- -------
         44     699084                 5.0 .. 37.18657261586156 29.277545928955078 .. 62.018354415893555                 5.0 .. 37.18657261586156 29.277545928955078 .. 62.018354415893555   17.13877296447754 .. 49.60246351587756         0     0    -inf     inf
        111     699078   33.88330268859863 .. 65.59163284301758  135.11531257629395 .. 167.1435375213623     37.948677458954 .. 69.65700761337295  135.11531257629395 .. 167.1435375213623  84.49930763244629 .. 116.36758518218994         0     0    -inf     inf
        217     699086 139.54998016357422 .. 170.55170822143555  243.45840644836426 .. 274.8020210266113  143.1504450626825 .. 174.15217312054384  243.45840644836426 .. 274.8020210266113 191.50419330596924 .. 222.67686462402344         0     0    -inf     inf
        322     699091   248.21048545837402 .. 278.352352142334   345.8145122528076 .. 376.4430561065674   251.7797389884793 .. 281.9216056724392   345.8145122528076 .. 376.4430561065674   297.0124988555908 .. 327.3977041244507         0     0    -inf     inf
        457     699080  350.76593017578125 .. 379.9788398742676   514.4993572235107 .. 544.0507583618164  355.0206685112119 .. 384.23357820969824    513.3126446738094 .. 542.864045812115     432.632643699646 .. 462.014799118042         0     0    -inf     inf
        553     699083   519.3173904418945 .. 546.6892547607422    539.7280474466841 .. 568.807982957441   519.3173904418945 .. 546.6892547607422    539.7280474466841 .. 568.807982957441   529.5227189442893 .. 557.7486188590916         2     2    -inf     inf
        569     699087   540.5885060803316 .. 569.6632521213445   549.8258158402863 .. 578.8447174665537   540.5885060803316 .. 569.6632521213445   549.8258158402863 .. 578.8447174665537     545.207160960309 .. 574.253984793949         3     3    -inf     inf

In addition, if reducing :doc:`deimos` data and slit-mask design matching is performed
(see :ref:`deimos:Slit-mask design matching`), a second
`astropy.io.fits.BinTableHDU`_ is written to disk.

.. code-block:: python

    Table(hdu[2].data)

    TRACEID TRACESROW     TRACELPIX          TRACERPIX      SLITID       SLITLOPT           SLITROPT         SLITRA      SLITDEC        SLITLEN       SLITWID  SLITPA ALIGN OBJID     OBJRA        OBJDEC
    int64    int64        float64            float64       int64        float64            float64         float64      float64        float64       float64 float64 int16 int64    float64      float64
    ------- --------- ------------------ ------------------ ------ ------------------- ------------------ ------------ ----------- ------------------ ------- ------- ----- ------ ------------ -----------
          0      2048 31.335221646243745 56.151084899902344 699084 -21.955958625056724  69.84514502686432 358.68126787  42.3347698 10.682000160217285     1.0     0.0     0 733790  358.6812625   42.334675
          1      2048  59.86090278625488 161.35546112060547 699078     72.846339998139  175.5490274983032 358.69310956 42.33811871 11.979000091552734     1.0     0.0     0 733784 358.69310417 42.33803333
          2      2048  164.9946002960205  269.1837863922119 699086  178.11184062306785 283.34493633572674 358.67168891 42.34139707 12.293000221252441     1.0     0.0     0 733792 358.67168333 42.34143889
          3      2048  272.9307098388672 371.00686264038086 699091   286.1885530241663   384.710234625776 358.65265564 42.34461342 11.529000282287598     1.0     0.0     0 733797    358.65265 42.34467778
          4      2048 374.73284912109375  538.8372478485107 699080  387.93013336193303  552.8797949320767 358.68945391 42.34933141 19.270000457763672     1.0     0.0     0 733786    358.68945 42.34819167
          5      2048   541.822380065918   563.737134898054 699083   555.6281596720763  576.8618516576732 358.68468468 42.35262222                4.0     4.0     0.0     1 733789 358.68467917 42.35262222
          6      2048  564.5937119658921  573.7892043748255 699087   576.9618516576732   586.896893857561 358.66447218 42.35271112                4.0     4.0     0.0     1 733793 358.66446667 42.35271111

See `Current SlitTrace Data Model`_ for a description of the columns.


Current SlitTrace Data Model
============================

Internally, the MasterSlits object is held in :class:`pypeit.slittrace.SlitTraceSet`,
which is a :class:`pypeit.datamodel.DataContainer`.


Here is its datamodel, which is written as an `astropy.io.fits.BinTableHDU`_.
In addition, if reducing :doc:`deimos` data and slit-mask design matching is performed,
another  `astropy.io.fits.BinTableHDU`_ is written to a fits extension named `MASKDEF_DESIGNTAB`.

.. include:: include/datamodel_slittrace.rst