.. _instruments:

=============
Spectrographs
=============

Overview
========

Below we describe all of the spectrographs that may
be reduced by PypeIt.  We also provide any suggested
tips for customizing the PypeIt file.

======================  =========   =======================================
PypeIt Name             Telescope   Instrument
======================  =========   =======================================
gemini_gmos_north_ham   Gemini      :doc:`gemini_gmos`-N spectrometer; Hamamatsu detector (R400, B600); Used since Feb 2017
gemini_gmos_north_e2v   Gemini      :doc:`gemini_gmos`-N spectrometer; E2V detector
gemini_gmos_south_ham   Gemini      :doc:`gemini_gmos`-S spectrometer; Hamamatsu detector (R400, B600)
gemini_gnirs            Gemini      GNIRS spectrometer
gemini_flamingos        Gemini      Gemini FLAMINGOS spectrometer
keck_kcwi               Keck        :doc:`keck_kcwi` slit-based IFU (BM, BH2)
keck_lris_blue          Keck        :doc:`lris` spectrometer; blue camera
keck_lris_red           Keck        :doc:`lris` spectrometer; red camera
keck_lris_red_orig      Keck        :doc:`lris` spectrometer; red camera + original detector
keck_mosfire            Keck        MOSFIRE spectrometer; J and Y gratings tested
keck_nires              Keck        NIRES spectrometer
keck_nirspec_low        Keck        NIRSPEC spectrometer; low-dispersion
keck_deimos             Keck        :doc:`deimos` spectrometer (600ZD, 830G, 1200G)
lbt_luci1               LBT         LUCI-I spectrometer
lbt_luci2               LBT         LUCI-II spectrometer
lbt_mods1r              LBT         MODS-I red spectrometer
lbt_mods1b              LBT         MODS-I blue spectrometer
lbt_mods2r              LBT         MODS-II red spectrometer
lbt_mods2b              LBT         MODS-II blue spectrometer
magellan_fire           Magellan    FIRE spectrometer; Echelle mode
magellan_fire_long      Magellan    FIRE spectrometer; Longslit high throughput mode
magellan_mage           Magellan    :doc:`mage` spectrometer
mdm_osmos               MDM         OSMOS spectrometer
mmt_mmirs               MMT         MMIRS spectrometer
mmt_binospec            MMT         BINSOSPEC spectrometer
not_alfosc              NOT         ALFOSC spectrometer (grisms 4, 19)
shane_kast_blue         Lick 3m     Kast dual spectrometer; blue camera
shane_kast_red          Lick 3m     Kast dual spectrometer; red camera
shane_kast_red_ret      Lick 3m     Kast dual spectrometer; red reticon
tng_dolores             TNG         DOLORES (LRS) spectrograph; LR-R
vlt_fors2               VLT         FORS2 spectrometer; only a few gratings
vlt_xshooter_uvb        VLT         :doc:`xshooter` spectrometer; UVB camera
vlt_xshooter_vis        VLT         :doc:`xshooter` spectrometer; VIS camera
vlt_xshooter_nir        VLT         :doc:`xshooter` spectrometer; NIR camera
wht_isis_blue           WHT         ISIS spectrometer; blue camera?
p200_dbsp_blue          P200        DBSP spectrograph; blue camera
p200_dbsp_red           P200        DBSP spectrograph; red camera
p200_tspec              P200        TripleSpec spectrograph
======================  =========   =======================================


List of Spectrographs
=====================

Instrument docs with additional details for running
PypeIt.


.. toctree::
   :caption: Spectrographs
   :maxdepth: 1

   new_spectrograph
   gemini_gmos
   deimos
   lris
   mage
   xshooter
