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

Near-IR
+++++++

If this is a near-IR instrument, you may wish to turn
off calibration steps.  See Gemini_GNIRS for an example.


