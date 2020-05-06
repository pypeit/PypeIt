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

This v1.0 release of ``PypeIt`` is designed to be used by both advanced
spectroscopists with prior data reduction expertise and new
astronomers. It is highly configurable and designed to be applied to
any standard slit-imaging spectrograph, and can accomodate longslit, multislit, as well 
as cross-dispersed echelle spectra. It has already enabled
several scientific publications [@hsyu2018; @eilers2020; @wang2020].

After the creation of a custom input/configuration file, the pipeline
runs end-to-end to convert raw spectroscopic images into calibrated,
science-ready spectra. It also includes scripts to flux calibrate and
combine multiple exposures, as well as software for performing telluric
corrections. ``PypeIt`` produces a series of
calibration-related outputs and includes scripts and automatically
generated plots for quality assurance inspection. The final outputs
are FITS files with rigid well-documented data models that hold the two-dimensional
(includes spatial information) and one-dimensional spectral
extractions.

It is our plan to expand ``PypeIt`` to include the majority of
spectrographs on the largest ground-based optical and near-infrared
telescopes, ideally with help from the broader community. Join us [on
GitHub](https://github.com/pypeit/PypeIt). By using this software, we
politely request that you abide by our [code of
conduct](https://pypeit.readthedocs.io/en/latest/codeconduct.html).

# Acknowledgements

We acknowledge intellectual contributions from TO BE FILLED IN. 

PypeIt has been financially supported by the University of California
Observatories. J.~F.~H. also acknowledges support from 
the University of California, Santa Barbara. During work on 
``PypeIt``,  R.~J.~C. was supported by a Royal Society University Research Fellowship, 
and acknowledges support from STFC (ST/P000541/1, ST/T000244/1).

# References
