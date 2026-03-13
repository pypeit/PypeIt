"""
Class that stores the file browsing element. Interacts with a backend interface
for file system inspection to allow for remote execution (if implemented).
"""

from __future__ import annotations

import os
import time
from typing import Dict, Iterable, List, Optional, Tuple

from ginga.misc import Bunch

from .backends import FileBrowserBackend
from .instruments import Instrument


class FileBrowserController:
    # Icon assets are injected after construction
    folderpb = None
    filepb = None
    fitspb = None

    def __init__(self, logger, settings, backend: FileBrowserBackend) -> None:
        """Initialize the FileBrowserController.

        Parameters
        ----------
        logger : logging.Logger
            Logger instance used for error and debug messages throughout the
            controller.
        settings : dict-like
            Configuration mapping providing browser settings such as
            ``"columns"`` (column definitions) and
            ``"max_rows_for_col_resize"`` (row count threshold above which
            column auto-resizing is skipped).
        backend : FileBrowserBackend
            Backend implementation responsible for filesystem operations
            (directory listing, stat calls, and FITS header extraction).
        """
        self.logger = logger
        self.settings = settings
        self.backend = backend

    def browse(
        self,
        path: str,
        instrument: Instrument,
        columns: Optional[List] = None,
        mode: str = "raw",
    ) -> Tuple[Dict[str, Bunch.Bunch], bool, str]:
        """Browse a directory and return a listing suitable for TreeView.

        Parameters
        ----------
        path : str
            Directory to list.
        instrument : Instrument
            Active instrument; used to read per-file metadata.
        columns : list of (str, str), optional
            Column definitions ``(display_name, attr_name)``.  When *None* the
            value from settings is used as a fallback.
        mode : str
            ``"raw"`` or ``"reduced"``; controls which instrument method is
            called to fetch FITS metadata.

        Returns
        -------
        listing : dict of {str : `~ginga.misc.Bunch.Bunch`}
            Tree dictionary keyed by filename suitable for use with a
            ``TreeView`` widget.  Each value is a ``Bunch`` containing file
            metadata fields and an ``icon`` pixbuf.
        resize : bool
            ``True`` when the number of rows is below the
            ``max_rows_for_col_resize`` settings threshold, indicating that
            column widths should be auto-resized after population.
        fullpath : str
            Glob-style path of the form ``<dirname>/*`` representing the
            directory that was listed.
        """
        if not os.path.isdir(path):
            raise ValueError(f"Invalid directory: {path}")

        if columns is None:
            columns = self.settings.get("columns")

        dirname = os.path.abspath(path)
        filelist = self.backend.list_dir(dirname)
        fullpath = os.path.join(dirname, "*")

        jumpinfo = list(
            map(lambda p: self._get_info(p, instrument, columns, mode), filelist)
        )
        listing, resize = self._makelisting(jumpinfo, columns)
        return listing, resize, fullpath

    def _get_info(
        self,
        path: str,
        instrument: Instrument,
        columns: List,
        mode: str = "raw",
    ) -> Bunch.Bunch:
        """Collect filesystem and instrument metadata for a single path.

        For plain FITS files the instrument backend is called to extract header
        metadata.  The same header extraction is performed for directories when
        *mode* is ``"reduced"``, allowing reduced-data directories to surface
        per-directory metadata in the listing.

        Parameters
        ----------
        path : str
            Absolute path to the file or directory to inspect.
        instrument : Instrument
            Active instrument; used to read per-file or per-directory header
            metadata via the backend.
        columns : list of (str, str)
            Column definitions as ``(display_name, attr_name)`` pairs.  Used
            to pre-populate missing attribute keys with the placeholder value
            ``"N/A"`` so that every returned ``Bunch`` has a consistent set of
            attributes regardless of whether metadata was available.
        mode : str, optional
            ``"raw"`` (default) or ``"reduced"``; controls which instrument
            method the backend calls to fetch FITS metadata.  In
            ``"reduced"`` mode, metadata is also attempted for directories.

        Returns
        -------
        bnch : `~ginga.misc.Bunch.Bunch`
            Bunch containing the following fields (at minimum):

            ``path`` : str
                Absolute path to the entry.
            ``name`` : str
                Basename of the entry.
            ``type`` : str
                One of ``"dir"``, ``"link"``, ``"fits"``, or ``"file"``.
            ``st_mode`` : int
                Raw permission/mode bits from ``os.stat``.
            ``st_mode_oct`` : str
                Octal string representation of *st_mode*.
            ``st_size`` : int
                File size in bytes.
            ``st_size_str`` : str
                String representation of *st_size*.
            ``st_mtime`` : float
                Modification time as a POSIX timestamp.
            ``st_mtime_str`` : str
                Human-readable modification time from :func:`time.ctime`.

            Additional keys are populated from the instrument header
            extraction when metadata is available, with column attribute names
            falling back to ``"N/A"`` when absent.
        """
        dirname, filename = os.path.split(path)
        name, ext = os.path.splitext(filename)
        ftype = "file"
        if os.path.isdir(path):
            ftype = "dir"
        elif os.path.islink(path):
            ftype = "link"
        elif filename.lower().endswith((".fits", ".fits.gz", ".fits.fz")):
            ftype = "fits"

        header_dict: Dict[str, object] = {}
        if ftype == "fits" or (ftype == "dir" and mode == "reduced"):
            try:
                header_dict = self.backend.get_header_info(path, instrument, mode=mode)
            except Exception as exc:
                self.logger.error(f"Error reading metadata for {path}: {exc}")
                header_dict = {}

        na_dict = {attrname: "N/A" for _colname, attrname in columns}
        bnch = Bunch.Bunch(na_dict)
        try:
            filestat = self.backend.stat_info(path)
            common = dict(
                path=path,
                name=filename,
                type=ftype,
                st_mode=filestat.st_mode,
                st_mode_oct=oct(filestat.st_mode),
                st_size=filestat.st_size,
                st_size_str=str(filestat.st_size),
                st_mtime=filestat.st_mtime,
                st_mtime_str=time.ctime(filestat.st_mtime),
            )
            if header_dict:
                common.update(header_dict)
            bnch.update(common)
        except OSError as exc:
            self.logger.error(f"Error getting file info for {path}: {exc}")
            bnch.update(dict(path=path, name=filename, type=ftype, st_mode=0, st_size=0, st_mtime=0))

        return bnch

    def _makelisting(
        self,
        jumpinfo: Iterable[Bunch.Bunch],
        columns: Optional[List] = None,
    ) -> Tuple[Dict[str, Bunch.Bunch], bool]:
        """Build the tree dictionary used to populate the file-browser widget.

        Iterates over the pre-collected file info bunches, attaches an
        appropriate icon pixbuf to each entry, and assembles the result into a
        dict keyed by filename.  Whether the caller should auto-resize table
        columns after population is determined by comparing the row count
        against the ``max_rows_for_col_resize`` settings value; auto-resize is
        skipped for large directories to avoid UI performance issues.

        Parameters
        ----------
        jumpinfo : iterable of `~ginga.misc.Bunch.Bunch`
            Sequence of file-info bunches as returned by :meth:`_get_info`.
        columns : list of (str, str), optional
            Column definitions ``(display_name, attr_name)``.  Currently
            accepted for interface consistency but not directly used in the
            listing construction.

        Returns
        -------
        tree_dict : dict of {str : `~ginga.misc.Bunch.Bunch`}
            Mapping from filename (``bnch.name``) to its info ``Bunch``.
            Each ``Bunch`` has an ``icon`` attribute set to the appropriate
            pixbuf for its file type.
        resize_flag : bool
            ``True`` when the number of entries is less than the
            ``max_rows_for_col_resize`` setting (default ``5000``), signalling
            that column auto-resizing should be performed.  ``False`` when the
            row count meets or exceeds the threshold.
        """
        def file_icon(bnch):
            """Return the pixbuf icon appropriate for the entry type.

            Parameters
            ----------
            bnch : `~ginga.misc.Bunch.Bunch`
                File-info bunch with a ``type`` attribute set to one of
                ``"dir"``, ``"fits"``, or any other string (treated as a
                generic file).

            Returns
            -------
            pixbuf : object
                ``FileBrowserController.folderpb`` for directories,
                ``FileBrowserController.fitspb`` for FITS files, or
                ``FileBrowserController.filepb`` for all other file types.
            """
            if bnch.type == "dir":
                return self.folderpb
            if bnch.type == "fits":
                return self.fitspb
            return self.filepb

        tree_dict: Dict[str, Bunch.Bunch] = {}
        for bnch in jumpinfo:
            icon = file_icon(bnch)
            bnch.setvals(icon=icon)
            entry_key = bnch.name
            if entry_key is None:
                raise Exception("No key for tuple")
            tree_dict[entry_key] = bnch

        n_rows = len(tree_dict)
        resize_table_columns = n_rows < self.settings.get("max_rows_for_col_resize", 5000)
        return tree_dict, resize_table_columns
