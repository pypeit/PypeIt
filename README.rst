
.. |PyPI| image:: https://img.shields.io/pypi/v/pypeit?label=PyPI&logo=pypi&logoColor=white
    :target: https://pypi.org/project/pypeit/#history

.. |pypi_downloads| image:: https://img.shields.io/pypi/dm/pypeit
    :target: https://pypi.org/project/pypeit/#files

.. |License| image:: https://img.shields.io/github/license/pypeit/PypeIt
   :target: https://github.com/pypeit/PypeIt/blob/release/LICENSE.rst

.. |CITests| image:: https://github.com/pypeit/PypeIt/workflows/CI%20Tests/badge.svg
    :target: https://github.com/pypeit/PypeIt/actions?query=workflow%3A"CI+Tests"

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

|docs| |CITests|

|DOI_latest| |JOSS| |arxiv|

PypeIt (pronounced "pipe it") is a Python package for semi-automated reduction
of astronomical spectroscopic data. Its algorithms build on decades-long
development of previous data reduction pipelines by the developers.

**For a complete description of PypeIt**, please see our `online documentation
<https://pypeit.readthedocs.io/en/stable/>`__.

----

.. _readme_usage:

Installation and Usage
++++++++++++++++++++++

Detailed installation instructions can be found `here
<https://pypeit.readthedocs.io/en/stable/installing.html>`__.  Briefly, after
creating a fresh python environment, simply type:

.. code-block:: bash

    pip install pypeit

PypeIt is designed to be used by both advanced spectroscopists with prior data
reduction expertise and astronomers with no prior experience of data reduction.
It is highly configurable and designed to be applied to any standard
slit-imaging spectrograph, including long-slit, multi-slit, as well as
cross-dispersed echelle spectra.

**The spectrographs that PypeIt supports are listed** `here
<https://pypeit.readthedocs.io/en/stable/spectrographs/spectrographs.html>`__.
Specifically, look `here
<https://pypeit.readthedocs.io/en/stable/spectrographs/spectrographs.html#instrument-specific-details>`__
for useful information about reducing data with certain instruments.

**Example data sets**: In addition to our primary code base, we maintain an
extensive `development suite
<https://github.com/pypeit/PypeIt-development-suite>`__ primarily used to
perform multiple layers of code testing, from basic unit tests to full
end-to-end tests of all our command-line scripts.  **If you are new to PypeIt**,
you are encouraged to learn how to use the code by finding and experimenting
with example data similar to your own in the ``RAW_DATA`` directory (organized
by instrument and configuration) of `this shared Google Drive folder
<https://drive.google.com/drive/folders/1oh19siB1-F0jjmY-F_jr73eA-TQYEiFW?usp=sharing>`__.

----

.. _community:

Community and Communication
+++++++++++++++++++++++++++

**Code of Conduct**: As a project, PypeIt is committed to fostering a
welcoming, diverse, and inclusive community.  As a member of this community you
are expected to read and follow our `Code of Conduct
<https://pypeit.readthedocs.io/en/stable/codeconduct.html>`__.

**Documentation**: We maintain extensive `online documentation
<https://pypeit.readthedocs.io/en/stable/>`__ that provides usage tutorials and
describes the PypeIt code, reduction and processing procedures, and output data
models.

**Real-time Communication**: We strongly encourage users to join 
our `PypeIt Users Slack <https://pypeit-users.slack.com>`__.  Developers are
available to answer questions and help troubleshoot issues.  All are welcome to
join using `this invitation link
<https://join.slack.com/t/pypeit-users/shared_invite/zt-3cderhn4g-XgFzv6mBWqxtVXKudc8W0w>`__;
please let us know (by submitting a GitHub issue) if the invitation link has
expired.  Please make sure to read the `guidelines
<https://docs.google.com/document/d/1Jo83RKstE8DaHxDCsGxsM3DC7sl153Hq9UGN3Piapuc/edit?usp=sharing>`__
before posting your question; the guidelines are also linked in the description
of the ``#guidelines`` channel.

**Bug Reports and Feature Requests**: If you encounter a bug, please first check
if this could be caused by user error; be sure to consult the online
documentation and post about your issue in the Users Slack (see above).  If the
problem persists (particularly if it is also experienced by other users),
please `submit a GitHub issue <https://github.com/pypeit/PypeIt/issues>`__.
When submitting the issue, please provide as much information as possible.  If
you haven't already, we may ask you to join the Users Slack for more efficient
communication about the issue.  Finally, you are also encouraged to `submit a
GitHub issue <https://github.com/pypeit/PypeIt/issues>`__ if you have a specific
feature request.

----

Contributing to PypeIt
++++++++++++++++++++++

.. NOTE: This is a copy of the CONTRIBUTING.rst file!  Try to keep them
   consistent.

We are excited to welcome your contributions to PypeIt!  We acknowledge 
contributions take many forms, including but not limited to participating in
discussions in our Users Slack Workspace; reporting issues to our GitHub
repository; submitting pull requests with small bug fixes, documentation
improvements, or large feature improvements; and participating in project
maintenance and governance.  All contributors are expected to follow our `Code
of Conduct <https://pypeit.readthedocs.io/en/stable/codeconduct.html>`__.

For direct contributions to the code, please see our `Development Guidelines
<https://pypeit.readthedocs.io/en/stable/dev/development.html>`__.  As mentioned
therein, communication between developers is key to ensuring efforts are
coordinated.  Before beginning any development activities, we would appreciate
communicating your intentions to the core development team, e.g., via the PypeIt
Users Slack Workspace.

For information regarding our governance structure and policies, please see the
`PypeIt Governance <https://pypeit.readthedocs.io/en/stable/governance.html>`__
documentation.

For a list of current contributors and project roles, please see our `PypeIt
Team <https://pypeit.readthedocs.io/en/stable/team.html>`__ listing.

----

Citation
++++++++

If you use PypeIt in your research, please cite the following two publications:

 - Prochaska et al. (2020, JOSS): `arXiv <https://ui.adsabs.harvard.edu/abs/2020arXiv200506505P/abstract>`__, `JOSS <https://joss.theoj.org/papers/10.21105/joss.02308>`__
 - Prochaska et al. (2020, Zenodo): `Zenodo <https://zenodo.org/records/17244617>`__

We provide the relevant `BibTeX entries
<https://pypeit.readthedocs.io/en/stable/index.html#pypeit-bibtex-entries>`__
for your convenience.  Note that the Zenodo BibTex entry is for version 1.0.0 of
the code; however, the Zenodo DOI is updated with each release of PypeIt.  When
you cite the Zenodo entry, you are encouraged to cite the specific PypeIt
version you used!  You are also encouraged to explicitly note the specific
PypeIt version used (e.g., 1.17.3) in the article text.

If there is no place to include the relevant citations in the
text of the publication, please include the following acknowledgement (provided
in latex and using the provided `BibTeX entries
<https://pypeit.readthedocs.io/en/stable/index.html#pypeit-bibtex-entries>`__):

.. code-block:: latex

    This research made use of \ttfamily{PypeIt} version
    1.17.3,\footnote{\url{https://pypeit.readthedocs.io/en/stable/}} a Python
    package for semi-automated reduction of astronomical slit-based spectroscopy
    \citep{pypeit:joss_pub, pypeit:zenodo}.

----

Funding
+++++++

PypeIt gratefully acknowledges funding from:

  * NASA ADAP (A20-0412, 20-1018)
  * NSF (TI-2346210, OAC-2410837)
  * JWST (JWST-AR-05464.001-A)
  * W\. M\. Keck Observatory
  * University of California Observatories

We also critically rely on important in-kind, open-source contributions from
the broader astronomical community.

