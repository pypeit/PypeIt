
.. |PyPI| image:: https://img.shields.io/pypi/v/pypeit?label=PyPI&logo=pypi&logoColor=white
    :target: https://pypi.org/project/pypeit/#history

.. |pypi_downloads| image:: https://img.shields.io/pypi/dm/pypeit
    :target: https://pypi.org/project/pypeit/#files

.. |License| image:: https://img.shields.io/github/license/pypeit/PypeIt
   :target: https://github.com/pypeit/PypeIt/blob/release/LICENSE.rst

.. |CITests| image:: https://github.com/pypeit/PypeIt/workflows/CI%20Tests/badge.svg
    :target: https://github.com/pypeit/PypeIt/actions?query=workflow%3A"CI+Tests"

.. |Coverage| image:: https://codecov.io/gh/PypeIt/pypeit/branch/release/graph/badge.svg
    :target: https://codecov.io/gh/PypeIt/pypeit

.. |docs| image:: https://readthedocs.org/projects/pypeit/badge/?version=latest
    :target: https://pypeit.readthedocs.io/en/latest/

.. |DOI_latest| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3743493.svg
   :target: https://doi.org/10.5281/zenodo.3743493

.. |JOSS| image:: https://joss.theoj.org/papers/10.21105/joss.02308/status.svg
   :target: https://doi.org/10.21105/joss.02308

.. |arxiv| image:: https://img.shields.io/badge/arxiv-2005.06505-black
   :target: https://arxiv.org/abs/2005.06505

.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org/

.. |forks| image:: https://img.shields.io/github/forks/pypeit/PypeIt?style=social
   :target: https://github.com/pypeit/PypeIt

.. |stars| image:: https://img.shields.io/github/stars/pypeit/PypeIt?style=social
   :target: https://github.com/pypeit/PypeIt

.. |github| image:: https://img.shields.io/badge/GitHub-PypeIt-brightgreen
   :target: https://github.com/pypeit/PypeIt


.. image:: https://raw.githubusercontent.com/pypeit/PypeIt/release/doc/_static/PypeIt_color_white_txt_black_background.png
    :target: https://github.com/pypeit/PypeIt
    :width: 485


PypeIt |forks| |stars|
======================

|github| |pypi| |pypi_downloads| |License|

|docs| |CITests| |Coverage| 

|DOI_latest| |JOSS| |arxiv|

PypeIt is a Python package for semi-automated reduction of astronomical
spectroscopic data. Its algorithms build on decades-long development of previous
data reduction pipelines by the developers. The reduction procedure - including
a complete list of the input parameters and available functionality - is
provided by our `online documentation
<https://pypeit.readthedocs.io/en/release/>`__.

PypeIt is designed to be used by both advanced spectroscopists with prior data
reduction expertise and astronomers with no prior experience of data reduction.
It is highly configurable and designed to be applied to any standard
slit-imaging spectrograph, including long-slit, multi-slit, as well as
cross-dispersed echelle spectra.  **The spectrographs that PypeIt can be used
with are listed** `here
<https://pypeit.readthedocs.io/en/release/spectrographs/spectrographs.html>`__.

Also note that `this link
<https://pypeit.readthedocs.io/en/release/spectrographs/spectrographs.html#instrument-specific-details>`__
contains some useful information about reducing data
with certain instruments that you might also find helpful.

In addition to our primary code base, we maintain an extensive `development
suite <https://github.com/pypeit/PypeIt-development-suite>`__ primarily used to
perform multiple layers of code testing, from basic unit tests to full
end-to-end tests of all our command-line scripts.  **If you are new to PypeIt,
you are encouraged to pull example data from the DevSuite for your instrument
when learning how to use the software.**

----

.. _community:

Community
+++++++++

As a project, PypeIt is committed to fostering a welcoming, diverse, and
inclusive community.  As a member of this community you are expected to read and
follow our `Code of Conduct
<https://pypeit.readthedocs.io/en/release/codeconduct.html>`__.

Along with our extensive `online documentation
<https://pypeit.readthedocs.io/en/release/>`__, we encourage the PypeIt user
base to communicate via our `PypeIt Users Slack <https://pypeit-users.slack.com>`__.
All are welcome to join using `this invitation link <https://join.slack.com/t/pypeit-users/shared_invite/zt-1kc4rxhsj-vKU1JnUA~8PZE~tPlu~aTg>`__.

If you find a bug (particularly one that is experienced by others in the Users
Slack) or have a feature request, please `submit an issue
<https://github.com/pypeit/PypeIt/issues>`__.

----

Citation
++++++++

If you use PypeIt in your research, please cite the following publications
(we provide the relevant `BibTeX entries
<https://pypeit.readthedocs.io/en/release/index.html#pypeit-bibtex-entries>`__
for your convenience):

 - Prochaska et al. (2020, JOSS): `arXiv <https://ui.adsabs.harvard.edu/abs/2020arXiv200506505P/abstract>`__, `JOSS <https://joss.theoj.org/papers/10.21105/joss.02308>`__
 - Prochaska et al. (2020, Zenodo): `Zenodo <https://ui.adsabs.harvard.edu/abs/2020zndo...3743493P/abstract>`__

If there is no place to include the relevant citations in the text of
the publication, please include the following acknowledgement
(provided in latex and using the provided `BibTeX entries
<https://pypeit.readthedocs.io/en/release/index.html#pypeit-bibtex-entries>`__):

.. code-block:: latex

    This research made use of \ttfamily{PypeIt},\footnote{\url{https://pypeit.readthedocs.io/en/latest/}}
    a Python package for semi-automated reduction of astronomical slit-based spectroscopy
    \citep{pypeit:joss_pub, pypeit:zenodo}.

----

Funding
+++++++

PypeIt receives direct funding from the following sources:

  * NASA ADAP (A20-0412, 20-1018)
  * W.M. Keck Observatory
  * University of California Observatories

We also rely on important in-kind contributions from individuals at
Caltech, the Multiple Mirror Observatory, and elsewhere.

----

Developers
++++++++++

PypeIt is an open-source, community developed package.  Astronomers are
encouraged to join the project and should review our `Code of Conduct
<https://pypeit.readthedocs.io/en/release/codeconduct.html>`__ and `Development
Guidelines <https://pypeit.readthedocs.io/en/release/dev/development.html>`__.
We would also appreciate if you contact the lead developers (JXP, JFH) before
beginning development activities.

The following persons have contributed substantially to the
development of PypeIt.

* J Xavier Prochaska (JXP)
* Joseph F. Hennawi (JFH)
* Kyle B. Westfall
* Ryan J. Cooke
* Feige Wang
* Tiffany Hsyu
* Frederick B. Davies
* Emanuele Paolo Farina
* Debora Pelliccia
* James Reichwein
* Milan Roberson
* Timothy Pickering
* Timothy Ellsworth-Bowers
* Gregory Simonian
* Heather Martin

