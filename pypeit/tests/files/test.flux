# User-defined fluxing parameters
[rdx]
   spectrograph = shane_kast_blue
[fluxcalib]
   std_file = spec1d_Feige66_KASTb_2015May20T041246.960.fits
   sensfunc = test_sensfunc.fits

flux read
  spec1d_J1217p3905_KASTb_2015May20T045733.560.fits test_flux.fits
flux end
