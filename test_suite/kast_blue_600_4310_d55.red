# This is a comment line

# Change the default settings
run ncpus 1
run spectrograph kast_blue
out verbose 2
out overwrite True
out sorted kast_blue_600_4310_d55

# Read in the data
data read
 /Users/xavier/PYPIT/Kast_blue/05192015/b*.fits.gz
data end

spect read
 fits calwin 12.
 pixflat number 3
 bias number 3
 trace number 3
 set standard b24.fits.gz 
spect end
