# PypeIt
[![PyPI](https://img.shields.io/pypi/v/pypeit?label=PyPI&logo=pypi&logoColor=white)](https://pypi.org/project/pypeit/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/pypeit?label=conda%20version)](https://anaconda.org/conda-forge/pypeit)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/pypeit?label=conda%20downloads)](https://anaconda.org/conda-forge/pypeit)

[![CI Tests](https://github.com/pypeit/PypeIt/workflows/CI%20Tests/badge.svg)](https://github.com/pypeit/PypeIt/actions?query=workflow%3A"CI+Tests")
[![Coverage (release)](https://codecov.io/gh/PypeIt/pypeit/branch/release/graph/badge.svg)](https://codecov.io/gh/PypeIt/pypeit)
[![Coverage (develop)](https://codecov.io/gh/PypeIt/pypeit/branch/develop/graph/badge.svg)](https://codecov.io/gh/PypeIt/pypeit)
[![Documentation Status](https://readthedocs.org/projects/pypeit/badge/?version=latest)](https://pypeit.readthedocs.io/en/latest/?badge=latest)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

The Python Spectroscopic Data Reduction Pipeline.  For
documentation visit:

http://pypeit.readthedocs.io

and/or see our HOWTO:

https://tinyurl.com/pypeit-howto

and/or join our PypeIt Users Slack
(the invite is recorded in this Issue:
https://github.com/pypeit/PypeIt/issues/676)

# Citation:

If you use ``PypeIt`` in your research, please cite the following
publications (BibTeX entries are provided below):

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02308/status.svg)](https://doi.org/10.21105/joss.02308)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3743493.svg)](https://doi.org/10.5281/zenodo.3743493)

If there is no place to include the relevant citations in the text of
the publication, please include the following acknowledgement
(provided in latex and using the provided BibTeX entries):

    This research made use of \ttfamily{PypeIt},\footnote{\url{https://pypeit.readthedocs.io/en/latest/}}
    a Python package for semi-automated reduction of astronomical slit-based spectroscopy
    \citep{pypeit:joss_pub, pypeit:zenodo}.

## BibTeX

    @ARTICLE{pypeit:joss_arXiv,
           author = {{Prochaska}, J. Xavier and {Hennawi}, Joseph F. and {Westfall}, Kyle B. and
             {Cooke}, Ryan J. and {Wang}, Feige and {Hsyu}, Tiffany and
             {Davies}, Frederick B. and {Farina}, Emanuele Paolo},
            title = "{PypeIt: The Python Spectroscopic Data Reduction Pipeline}",
          journal = {arXiv e-prints},
         keywords = {Astrophysics - Instrumentation and Methods for Astrophysics},
             year = 2020,
            month = may,
              eid = {arXiv:2005.06505},
            pages = {arXiv:2005.06505},
    archivePrefix = {arXiv},
           eprint = {2005.06505},
     primaryClass = {astro-ph.IM},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2020arXiv200506505P},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

    @article{pypeit:joss_pub,
        doi = {10.21105/joss.02308},
        url = {https://doi.org/10.21105/joss.02308},
        year = {2020},
        publisher = {The Open Journal},
        volume = {5},
        number = {56},
        pages = {2308},
        author = {J. Xavier Prochaska and Joseph F. Hennawi and Kyle B. Westfall and Ryan J. Cooke and Feige Wang and Tiffany Hsyu and Frederick B. Davies and Emanuele Paolo Farina and Debora Pelliccia},
        title = {PypeIt: The Python Spectroscopic Data Reduction Pipeline},
        journal = {Journal of Open Source Software}
    }

    @MISC{pypeit:zenodo,
           author = {{Prochaska}, J. Xavier and {Hennawi}, Joseph and {Cooke}, Ryan and
             {Westfall}, Kyle and {Wang}, Feige and {EmAstro} and {Tiffanyhsyu} and
             {Wasserman}, Asher and {Villaume}, Alexa and {Marijana777} and
             {Schindler}, JT and {Young}, David and {Simha}, Sunil and
             {Wilde}, Matt and {Tejos}, Nicolas and {Isbell}, Jacob and
             {Fl{\"o}rs}, Andreas and {Sandford}, Nathan and {Vasovi{\'c}}, Zlatan and
             {Betts}, Edward and {Holden}, Brad},
            title = "{pypeit/PypeIt: Release 1.0.0}",
             year = 2020,
            month = apr,
              eid = {10.5281/zenodo.3743493},
              doi = {10.5281/zenodo.3743493},
          version = {v1.0.0},
        publisher = {Zenodo},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2020zndo...3743493P},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

# Contribute

We encourage anyone to help us develop the `PypeIt` code base to better
suit your needs and to improve its algorithms. If you do so, please
follow our [Development
Guidlines](https://pypeit.readthedocs.io/en/latest/development.html)

In particular, please note our [Code of
Conduct](https://pypeit.readthedocs.io/en/latest/codeconduct.html).


# Instruments Served
* Bok/B&C
* Gemini/GNIRS
* Gemini/GMOS
* Gemini/FLAMINGOS 2
* GTC/OSIRIS
* Lick/Kast
* Magellan/MagE
* Magellan/Fire
* MMT/BinoSpec (270 and 600 tested)
* MMT/MMIRS (HK_zJ, J_zJ, and K_K tested)
* MMT/Blue Channel (300 tested)
* MDM/OSMOS
* Keck/DEIMOS (600ZD, 830G, 1200G)
* Keck/KCWI (BM, BH2)
* Keck/LRIS
* Keck/MOSFIRE  (J and Y gratings tested)
* Keck/NIRES
* Keck/NIRSPEC (low-dispersion)
* LBT/Luci-I, Luci-II
* LBT/MODS (beta)
* LDT/DeVeny
* Lick/APF (planned)
* NOT/ALFOSC (grism4)
* VLT/X-Shooter
* VLT/FORS2  (300I, 300V)
* WHT/ISIS
* P200/DBSP (316/7500 on red arm, 600/4000 on blue arm)
* P200/TripleSpec

# Requirements

(see `setup.cfg` or `environment.yml`)

* python
* numpy
* scipy
* matplotlib
* astropy
* ginga
* h5py
* future
* PyYAML
* linetools
* IPython
* scikit-learn
* configobj


# License (BSD-3)

(see `LICENSE.rst`)

Copyright (c) 2018-2019, PypeIt Developers All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

 - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

 - Neither the name of the Astropy Team nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
