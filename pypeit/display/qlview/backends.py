from __future__ import annotations

import glob
import os
from typing import Dict, List, Protocol, TYPE_CHECKING

if TYPE_CHECKING:
    from .instruments import Instrument

from pypeit.scripts import ql


class FileBrowserBackend(Protocol):
    def list_dir(self, path: str) -> List[str]:
        ...

    def stat_info(self, path: str) -> os.stat_result:
        ...

    def get_header_info(self, path: str, instrument: "Instrument", mode: str = "raw") -> Dict[str, object]:
        ...


class ReductionBackend(Protocol):
    def submit(self, args: List[str]) -> None:
        ...


class LocalFileBrowserBackend:
    def list_dir(self, path: str) -> List[str]:
        fullpath = os.path.join(os.path.abspath(path), "*")
        filelist = list(glob.glob(fullpath))
        filelist.sort(key=lambda s: s.lower())
        filelist.insert(0, os.path.join(os.path.abspath(path), ".."))
        return filelist

    def stat_info(self, path: str) -> os.stat_result:
        return os.stat(path)

    def get_header_info(self, path: str, instrument: "Instrument", mode: str = "raw") -> Dict[str, object]:
        if mode == "reduced":
            return instrument.get_reduced_info(path)
        return instrument.get_raw_info(path)


class RemoteFileBrowserBackend:
    def __init__(self, base_url: str) -> None:
        self.base_url = base_url

    def list_dir(self, path: str) -> List[str]:
        raise NotImplementedError("Remote backend not implemented")

    def stat_info(self, path: str) -> os.stat_result:
        raise NotImplementedError("Remote backend not implemented")

    def get_header_info(self, path: str, instrument: "Instrument", mode: str = "raw") -> Dict[str, object]:
        raise NotImplementedError("Remote backend not implemented")


class LocalReductionBackend:
    def submit(self, args: List[str]) -> None:
        parsed = ql.QL.get_parser().parse_args(args)
        ql.QL.main(parsed)


class RemoteReductionBackend:
    def __init__(self, base_url: str) -> None:
        self.base_url = base_url

    def submit(self, args: List[str]) -> None:
        raise NotImplementedError("Remote backend not implemented")
