---
title: 'PypeIt: The Python Spectroscopic Data Reduction Pipeline'
tags:
  - Python
  - astronomy
  - data reduction
  - spectroscopy
authors:
  - name: J. Xavier Prochaska
    orcid: 0000-0002-7738-6875
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Joseph F. Hennawi
    orcid: 0000-0002-7054-4332
    affiliation: 3
  - name: Kyle B. Westfall
    orcid: 0000-0003-1809-6920
    affiliation: 4
  - name: Ryan J. Cooke
    orcid: 0000-0001-7653-5827
    affiliation: 5
  - name: Feige Wang
    orcid: 0000-0002-7633-431X
    affiliation: 6
  - name: Tiffany Hsyu
    orcid: 0000-0002-0462-3139
    affiliation: 1
  - name: Emanuele Paolo Farina
    orcid: 0000-0002-6822-2254
    affiliation: 7
affiliations:
 - name: University of California, Santa Cruz
   index: 1
 - name: Kavli Institute for the Physics and Mathematics of the Universe
   index: 2
 - name: University of California, Santa Barbara
   index: 3
 - name: University of California Observatories
   index: 4
 - name: Durham University, UK
   index: 5
 - name: Steward Observatory, University of Arizona
   index: 6
 - name: Max Planck Institut f\"ur Astrophysik
   index: 7
date: 08 April 2020
bibliography: paper.bib
---

# Summary

``PypeIt`` is a Python package for semi-automated reduction of
astronomical, spectroscopic data. Its algorithms build on
decades-long development of previous data reduction pipelines by the
developers [@mike; @mase]. The reduction procedure - including a
complete list of the input parameters and available functionality -
is provided as online documentation hosted by [Read the
Docs](https://pypeit.readthedocs.io), which is regularly updated.
In what follows, we provide a brief description of the algorithms,
but refer the interested reader to the online documentation. 

Release v1.0.3 serves the following spectrographs:
Gemini/GNIRS, Gemini/GMOS, Gemini/FLAMINGOS 2, Lick/Kast, Magellan/MagE,
Magellan/Fire, MDM/OSMOS, Keck/DEIMOS (600ZD, 830G, 1200G), Keck/LRIS,
Keck/MOSFIRE (J and Y gratings tested), Keck/NIRES, Keck/NIRSPEC
(low-dispersion), LBT/Luci-I, Luci-II, LBT/MODS (beta), NOT/ALFOSC (grism4),
VLT/X-Shooter (VIS, NIR), VLT/FORS2 (300I, 300V), WHT/ISIS.

This v1.0 release of ``PypeIt`` is designed to be used by both advanced
spectroscopists with prior data reduction expertise and astronomers with
no prior experience of data reduction. It is highly configurable and
designed to be applied to any standard slit-imaging spectrograph, and
can accomodate longslit, multislit, as well as cross-dispersed echelle
spectra. It has already enabled several scientific publications
[@hsyu2018; @eilers2020; @wang2020].

After the creation of a custom input/configuration file, the pipeline
runs end-to-end to convert raw spectroscopic images into calibrated,
science-ready spectra. In what follows, we describe several key steps
of the data reduction procedure:

(1) bias, dark, bpm,

(2) slits traces

(3) Master arc, wavelength calibration, identify GUI

(4) tilts image

(5) pixelflat, spatial slit map

(6) spatial and spectral flexure correction

(7) reference frame heliocentric, barycentric

(8) advanced sky subtraction routines, including a modified version of the Kelson algorithm.

(9) boxcar + optimal extraction

PypeIt also includes scripts to flux calibrate and
combine multiple exposures, as well as software for performing telluric
corrections. ``PypeIt`` produces a series of
calibration-related outputs and includes scripts and automatically
generated plots for quality assurance inspection. The final outputs
are FITS files with rigid well-documented data models that hold the two-dimensional
(includes spatial information) and one-dimensional spectral
extractions.

It is our plan to expand ``PypeIt`` to include the majority of
spectrographs on the largest ground-based optical and near-infrared
telescopes, ideally with help from the broader community. We are
currently working towards implementing the following additional
spectrographs:
Keck/DEIMOS (all gratings)
Keck/KCWI,
Keck/MOSFIRE (all gratings),
Keck/NIRSPEC (new detector + high resolution),
Magellan/IMACS,
MMT/BinoSpec. We are also open to receiving
requests to support additional spectroscopic
instrumentation.
Join us [on GitHub](https://github.com/pypeit/PypeIt). By using this software, we
politely request that you abide by our [code of
conduct](https://pypeit.readthedocs.io/en/latest/codeconduct.html).

# Acknowledgements

We acknowledge intellectual contributions from Scott Burles,
Rob Simcoe, and David Schlegel.

PypeIt has been financially supported by the University of California
Observatories. J.~F.~H. also acknowledges support from 
the University of California, Santa Barbara. During work on 
``PypeIt``,  R.~J.~C. was supported by a Royal Society University Research Fellowship, 
and acknowledges support from STFC (ST/P000541/1, ST/T000244/1).

# References
