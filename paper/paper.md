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
    affiliation: 3
affiliations:
 - name: University of California, Santa Cruz
   index: 1
 - name: Kavli Institute for the Physics and Mathematics of the Universe
   index: 2
 - name: University of California, Santa Barbara
   index: 3
date: 08 April 2020
bibliography: paper.bib

# Summary

`PypeIt` is a Python package for semi-automated reduction of 
astronomical, spectroscopic data. 
Its algorithms build on decades-long development of previous
data reduction pipelines by the developers [@mike; @mase].

After the creation of a custom input/configuration file,
the pipeline runs end-to-end to convert raw spectroscopic images
into calibrated, science-ready spectra.
It also includes scripts to flux and combine multiple exposures.
`PypeIt` produces a series of calibration-related outputs and includes
scripts for quality assurance inspection.  The final outputs
are FITS files with rigid data models that hold the
two-dimensional (includes spatial information) and
one-dimensional spectral extractions.

Release v1.0 of `PypeIt` is designed to be used by both advanced 
spectroscopists with prior data-reduction expertise and new astronomers. 
It is highly configurable and designed to be applied to any 
standard spectrograph.
It has already enabled several scientific publications 
[@hsyu2018; @eilers2020]. 

It is our plan to expand PypeIt to include the majority of spectrographs
on the largest ground-based telescopes, ideally with help from
the broader community.


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"


# Acknowledgements

We acknowledge contributions from TO BE FILLED IN. 

PypeIt has been financially supported by 
the University of California Observatories.

# References