.. highlight:: rest

=================
PyPI Instructions
=================

Setup
-----

 - You will need to install `twine` via pip on your machine.
 - Generate a .pypirc file in your home directory that looks like this::

    [distutils]
    index-servers = pypi

    [pypi]
    repository = https://upload.pypi.org/legacy/
    username = pypeit
    password = [ask for this]



Do it
-----

Am following instructions for `twine <https://pypi.org/project/twine>`_

Step 1 (Build)::

    python setup.py sdist bdist_wheel

Step 2 (Test)::

    twine upload --repository-url https://test.pypi.org/legacy/ dist/*

Step 3 (For real)::

    twine upload dist/*
