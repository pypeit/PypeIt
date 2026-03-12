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
        dirname, filename = os.path.split(path)
        name, ext = os.path.splitext(filename)
        ftype = "file"
        if os.path.isdir(path):
            ftype = "dir"
        elif os.path.islink(path):
            ftype = "link"
        elif ext.lower() == ".fits":
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
        def file_icon(bnch):
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
