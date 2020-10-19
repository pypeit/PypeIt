.. include:: ../include/links.rst

.. _conda_forge:

How to manually update the PypeIt conda-forge feedstock (what to do if the bots are broken)
===========================================================================================

1. Fork `conda-forge/pypeit-feedstock <https://github.com/conda-forge/pypeit-feedstock>`_.
2. Install required tools

.. code-block:: bash

    conda install -c conda-forge conda-smithy

3. Manually update the ``version`` variable in recipe/meta.yaml
4. Update the ``source.sha256`` field in recipe/meta.yaml with the sha256 hash available on
   `pypi <https://pypi.org/project/pypeit/#files>`_

.. code-block:: yaml

    source:
        sha256: new_hash_here

5. Reset the ``build.number`` field to 0 in recipe/meta.yaml

.. code-block:: yaml

    build:
        number: 0

6. Commit your changes
7. Rerender the feedstock

.. code-block:: bash

    conda smithy rerender -c auto

8. Update the dependencies in recipe/meta.yaml if PypeIt's dependencies have changed.
9. Open a PR with these changes.
10. Merge the PR after all tests pass.
