
.. include:: ../include/links.rst

.. _new_script:

*****************************
Developing New PypeIt Scripts
*****************************

All of the PypeIt executable scripts are located in the ``pypeit/scripts``
directory, and they all have roughly the same structure:

.. code-block:: python

    from pypeit.scripts import scriptbase

    class NewScript(scriptbase.ScriptBase):

        @classmethod
        def get_parser(cls, width=None):
            parser = super().get_parser(description='A new PypeIt script', width=width)
            ...

        @staticmethod
        def main(args):
            ...

The important components of the scripts are:

 * To ease installation and documentation of the scripts, they all must use
   :class:`~pypeit.scripts.scriptbase.ScriptBase` as their base class.  Note
   that the class has no instantiation method.

 * The :func:`~pypeit.scripts.scriptbase.ScriptBase.get_parser` function is used
   to return an instance of `argparse.ArgumentParser`_ used to parse the
   command-line arguments.

 * The :func:`~pypeit.scripts.scriptbase.ScriptBase.main` function performs the
   primary operations of the script.  It needs no return value, but a return is
   not prohibited (see, e.g.,
   :class:`~pypeit.scripts.chk_for_calibs.ChkForCalibs`).

 * Each file in the ``pypeit/scripts`` directory should only contain *one*
   script class.  The :func:`~pypeit.scripts.scriptbase.ScriptBase.name`
   function sets the name of the script to ``pypeit_{module}`` by default, where
   ``{module}`` is the name of the file.  E.g., if the file name of the new
   script is ``new_script.py`` the executable installed will be
   ``pypeit_new_script``.  This can be altered by overriding the base class
   ``name`` function (see :class:`~pypeit.scripts.run_pypeit.RunPypeIt`).

The base class, :class:`~pypeit.scripts.scriptbase.ScriptBase`, provides the
common entry point function
(:func:`~pypeit.scripts.scriptbase.ScriptBase.entry_point`) used during
installation.  To ensure that the script is properly installed by `pip`_, you
need to add this entry point to the ``setup.cfg`` file.  All of the PypeIt
scripts are listed in the ``[options.entry_points]`` group.  To add your script,
you enter a new line with the following format:

.. code-block:: ini

    pypeit_new_script = pypeit.scripts.new_script:NewScript.entry_point

Lastly, you should add the import of the script module to the
``pypeit/scripts/__init__.py`` file; i.e., add:

.. code-block:: python

    from pypeit.scripts import new_script

Note that the script files in the ``pypeit/scripts`` directory:

 * Should not be executable (i.e., no xs in their permissions)
 * Should not start with an env statement; i.e., ``#!/usr/bin/env python``
 * Should not end with the ``if __name__ == '__main__':`` block

Creating the executables from the raw script files is all handled by `pip`_
installing PypeIt.  To ensure the script is installed, from the top-level
directory run, e.g.:

.. code-block:: console

    pip install -e ".[dev]"

If the new script doesn't appear in your path after running this, you may need
to uninstall (``pip uninstall pypeit``) and reinstall using the command above.


