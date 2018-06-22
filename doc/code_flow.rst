.. highlight:: rest

=========
Code Flow
=========

This describes the standard code flow of PYPIT

ARMS
====

Multi-slit and longslit reductions.

===============  ============= ============= ============================================ ===========================
Step             Class         Internals     Outputs                                      QA
===============  ============= ============= ============================================ ===========================
Setup            PypitSetup    fitstbl       keck_lris_red_setup_A.fits
..                             setup_dict    setup_files/keck_lris_red_2018-Jun-19.setups
..                                           keck_lris_red_setup_A.pypit
Bias             BiasFrame     msbias        MasterBias_A_02_aa.fits
ArcImg           ArcImage      msarc         MasterArc_A_02_aa.fits
Bad Pixel Mask   BPMImage      msbpm
Pixel location   ---           pixlocn
Trace Slits      TraceSlit     tslits_dict   MasterTrace_A_02_aa.fits.gz                  Slit_Trace_A_02_aa.png
..                                           MasterTrace_A_02_aa.json
1D Wave Calib    WaveCalib     wv_calib      MasterWaveCalib_A_02_aa.json                 Arc_1dfit_A_02_aa_S0000.png
Wave Tilts       WaveTilts     mstilts       MasterTilts_A_02_aa.fits                     Arc_tilts_A_02_aa_S0000.png
Pixel flat       FlatField     mspixflatnrm  MasterFlatField_A_02_aa.fits
..
Process Science  ScienceImage  sciframe      spec2d_basename.fits
Global skysub    ScienceImage  global_sky
Find objects     ScienceImage  tracelist
Extraction       ScienceImage  specobjs
Flexure          arflex        flex_list                                                  basename_flex_sky.png
..                                                                                        basename_flex_corr.png
===============  ============= ============= ============================================ ===========================
