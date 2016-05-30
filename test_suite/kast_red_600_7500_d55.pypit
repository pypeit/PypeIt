# This is a comment line

# Change the default settings
run ncpus 1
run spectrograph kast_red
out verbose 2
out overwrite True
out sorted kast_red_600_7500_d55

# Reduce
#arc calibrate id_pix 220.63808792,417.17122925,562.95948074,669.07189712,884.38191292
#arc calibrate id_wave 5868.4165,6334.4276,6678.2766,6929.4672,7438.8981
trace orders pcatilt 1,1,1
trace orders tilts spca

# Read in the data
data read
 ~/Dropbox/PYPIT/TEST_SUITES/Kast_red/600_7500_d55/r1*.fits.gz
data end

spect read
 fits calwin 12.
 trace number 3
 pixflat number 3
 bias number 3
 set standard r136.fits.gz
 set arc r124.fits.gz
spect end
