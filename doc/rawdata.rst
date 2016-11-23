.. highlight:: rest

********
Raw Data
********

Raw Files
=========

The raw files may be located anywhere on the harddrive.
We recommend you organize them by date of observation.

You will need to provide the full path to the files
and a filename root (unless you wish to specify them
one by one!).  Here is an example path+root::

    /Users/xavier/Keck/LRIS/data/2016apr06/Raw/LB.20160406.

The code will work on all files with this root.

Compression
===========

The files can be gzipped with a .gz extension
but we warn that Python's FITS reader works
considerably slower on compressed files.

Therefore, you may wish to uncompress while reducing
and compress again afterwards.


