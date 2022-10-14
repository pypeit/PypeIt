.. include:: ../README.rst

----

.. _bibtex:

PypeIt BibTeX Entries
+++++++++++++++++++++

.. code-block:: latex

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

.. toctree::
   :caption: Welcome!
   :maxdepth: 1
   :hidden:

   codeconduct
   installing
   whatsnew
   known_failure_modes

.. toctree::
   :caption: Quickstart & Tutorials
   :maxdepth: 1
   :hidden:

   cookbook
   tutorials/tutorials
   quicklook
   Example Data <dev_suite>
   reduction_tips

.. toctree::
   :caption: Core Processing
   :maxdepth: 1
   :hidden:

   setup
   pypeit_par
   Execution <running>
   outputs
   QA <qa>

.. toctree::
   :caption: Processing Details
   :maxdepth: 1
   :hidden:

   calibrations/calibrations
   2d_combine
   object_finding
   skysub
   extraction
   A-B_differencing

.. toctree::
   :caption: Further Processing
   :maxdepth: 1
   :hidden:

   fluxing
   coadd1d
   telluric
   collate1d
   coadd2d
   coadd3d

.. toctree::
   :caption: Reference
   :maxdepth: 1
   :hidden:

   pypeit_par
   input_files
   pypeit_file
   conventions
   spectrographs/spectrographs
   detectors
   mosaic
   frametype
   scripts
   bitmasks

.. toctree::
   :caption: For Developers
   :maxdepth: 1
   :hidden:

   Development Guidelines <dev/development>
   Adding a New Spectrograph <dev/new_spectrograph>
   API <api/modules>

