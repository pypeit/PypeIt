.. highlight:: rest

****************
New Spectrograph
****************

Here are notes on how to build a new spectrograph
from scratch or to add a new mode.

Entirely New
============


#.  Build a new name_of_spectrograph.py file in pypeit.spectrograph
#.  Fuss with the Detector object; one per detector
    - Set datasec, oscansec in the *raw* frame, i.e. as viewed on Ginga
    - Or generate a custom reader if these are variable
#.  Set custom parameters


