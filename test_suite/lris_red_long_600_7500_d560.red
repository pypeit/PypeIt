# This is a comment line

# Change the default settings
run ncpus 1
run spectrograph lris_red
out verbose 2
out overwrite True
out sorted lris_r_ls_600_7500_d560

# Reduce
reduce usebias bias      # How to subtract the detector bias (bias, overscan, dark, none), you can also specify a master calibrations file if it exists.
arc calibrate id_pix 140.93677186,542.27337912,739.04374828,1152.38285967,1767.86702435
arc calibrate id_wave 5771.210,6404.018,6718.8974,7386.014,8379.9093
trace disp direction 0   # Manually specify the dispersion direction (0 for row, 1 for column) 0 for column, 1 for row????
trace orders pcatilt 1,1,1
trace orders tilts spca
reduce bgsubtraction method bspline
pixflat comb method median
pixflat comb rej_level [10.0,10.0]
pixflat norm recnorm False

# Read in the data
data read
 /Users/xavier/Dropbox/PYPIT/LRIS_red_xmps_test/LR*.fits
data end

spect read
 fits calwin 12.
 pixflat number 5
 bias number 3
 arc number 1
 trace number 5
 set arc LR.20160216.05529.fits
 set arc LR.20160216.05589.fits
 set arc LR.20160216.05649.fits
 set arc LR.20160216.05709.fits
 set bias LR.20160216.07348.fits
 set bias LR.20160216.07412.fits
 set bias LR.20160216.07470.fits
 set bias LR.20160216.07529.fits
 set bias LR.20160216.07587.fits
 set bias LR.20160216.07646.fits
 set bias LR.20160216.07703.fits
 set bias LR.20160216.07762.fits
 set bias LR.20160216.07820.fits
 set bias LR.20160216.07878.fits
 set bias LR.20160216.07937.fits
 set trace LR.20160216.13991.fits
 set trace LR.20160216.14090.fits
 set trace LR.20160216.14167.fits
 set trace LR.20160216.14244.fits
 set trace LR.20160216.14322.fits
 set trace LR.20160216.14399.fits
 set pixflat LR.20160216.13991.fits
 set pixflat LR.20160216.14090.fits
 set pixflat LR.20160216.14167.fits
 set pixflat LR.20160216.14244.fits
 set pixflat LR.20160216.14322.fits
 set pixflat LR.20160216.14399.fits
 set standard LR.20160216.17613.fits
 set science LR.20160216.40478.fits
spect end
