"""
Command handlers run by the :class:`~pypeit.display.pypeit_remote.PypeItRemote`
ginga plugin.

Each handler runs one PypeIt display command inside the ginga process.  It is
called on a ginga worker thread as ``handler(fv, args, job)``, where ``fv`` is
the ginga reference viewer shell, ``args`` is the `argparse.Namespace`_ parsed
by the script that sent the command, and ``job`` is the
:class:`~pypeit.display.remote_jobs.RemoteJob` running it.

A handler must:

    - do file I/O and computation directly, on the worker thread;
    - make every change to the viewer through ``fv.gui_call``;
    - call ``job.check_cancel()`` between steps so that a newer command can
      supersede it;
    - raise :class:`~pypeit.PypeItError` for errors the user should see in the
      terminal.

Warnings logged with the PypeIt logger by the worker thread are reported to
the script that sent the command.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import PypeItError
from pypeit import specobjs
from pypeit.display import ginga_draw


def get_handlers():
    """
    Return the command handlers known to the ``PypeItRemote`` plugin.

    The command names match the names used by the scripts that send them
    (e.g., ``'show_2dspec'`` for ``pypeit_show_2dspec``).

    Returns
    -------
    :obj:`dict`
        Handler functions keyed by command name.
    """
    return {'show_1dspec': show_1dspec}


def clear_workspace(fv):
    """
    Delete all channels in the viewer.

    Must be called on the ginga GUI thread.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    """
    ginga_draw.clear_channels(fv, list(fv.get_channel_names()))


def spec1d_extension(sobjs, exten=1, obj=None, extract='OPT'):
    """
    Select the spectrum to show from a spec1d file.

    Parameters
    ----------
    sobjs : :class:`~pypeit.specobjs.SpecObjs`
        Spectra read from the spec1d file.
    exten : :obj:`int`, optional
        1-based number of the spectrum to show.  Ignored if ``obj`` is
        provided.
    obj : :obj:`str`, optional
        Name of the spectrum to show, e.g. ``SPAT0424-SLIT0000-DET01``.
    extract : :obj:`str`, optional
        Extraction to show, ``'OPT'`` or ``'BOX'``.

    Returns
    -------
    :obj:`int`
        0-based index of the spectrum in ``sobjs``.

    Raises
    ------
    :class:`~pypeit.PypeItError`
        Raised if the spectrum does not exist or does not have the requested
        extraction.
    """
    nspec = len(sobjs)
    if obj is None:
        indx = exten - 1
        if indx < 0 or indx >= nspec:
            raise PypeItError(f'Extension {exten} is not available; the file has {nspec} spectra.')
    else:
        match = np.where(sobjs.NAME == obj)[0]
        if len(match) == 0:
            raise PypeItError(f'No spectrum named {obj}; use --list to see the names.')
        indx = int(match[0])
    if extract == 'OPT' and sobjs[indx]['OPT_WAVE'] is None:
        raise PypeItError('Spectrum not extracted with OPT.  Try --extract BOX')
    return indx


def _start_spec1dview(fv, chname, params):
    """
    Clear the viewer and start the Spec1dView plugin in a new channel.

    Must be called on the ginga GUI thread.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    chname : :obj:`str`
        Name of the channel to create.
    params : :obj:`dict`
        Keyword arguments for ``Spec1dView.set_params``.
    """
    clear_workspace(fv)
    fv.add_channel(chname)
    fv.call_local_plugin_method(chname, 'Spec1dView', 'set_params', [], params)
    fv.start_local_plugin(chname, 'Spec1dView')
    fv.change_channel(chname)


def show_1dspec(fv, args, job):
    """
    Show a spectrum from a spec1d file with the Spec1dView plugin.

    Handler for ``pypeit_show_1dspec``; see
    :class:`~pypeit.scripts.show_1dspec.Show1DSpec` for the arguments.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    args : `argparse.Namespace`_
        Arguments parsed by ``pypeit_show_1dspec``.
    job : :class:`~pypeit.display.remote_jobs.RemoteJob`
        The job running the command.
    """
    sobjs = specobjs.SpecObjs.from_fitsfile(args.file, chk_version=False)
    indx = spec1d_extension(sobjs, exten=args.exten, obj=args.obj, extract=args.extract)
    job.check_cancel()
    chname = 'Spec1d'
    params = dict(ext=indx, extraction=args.extract, masked=args.masked, fluxed=args.flux)
    fv.gui_call(_start_spec1dview, fv, chname, params)
    job.check_cancel()
    fv.load_file(args.file, chname=chname, wait=True)
