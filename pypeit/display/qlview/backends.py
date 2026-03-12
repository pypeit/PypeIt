from __future__ import annotations

import glob
import os
import types
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
    """File-browser backend that delegates to a remote :mod:`HTTPserver`."""

    def __init__(self, base_url: str) -> None:
        self.base_url = base_url.rstrip("/")

    def list_dir(self, path: str) -> List[str]:
        import requests
        resp = requests.get(
            f"{self.base_url}/api/list_dir",
            params={"path": path},
            timeout=30,
        )
        resp.raise_for_status()
        return resp.json()["files"]

    def stat_info(self, path: str) -> os.stat_result:
        import requests
        resp = requests.get(
            f"{self.base_url}/api/stat_info",
            params={"path": path},
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
        # os.stat_result cannot be constructed directly; return a namespace that
        # exposes the fields consumed by FileBrowserController.
        return types.SimpleNamespace(
            st_mode=data["st_mode"],
            st_size=data["st_size"],
            st_mtime=data["st_mtime"],
            st_uid=data.get("st_uid", 0),
            st_gid=data.get("st_gid", 0),
        )

    def get_header_info(self, path: str, instrument: "Instrument", mode: str = "raw") -> Dict[str, object]:
        import requests
        resp = requests.get(
            f"{self.base_url}/api/header_info",
            params={
                "path": path,
                "instrument": instrument.__class__.__name__,
                "mode": mode,
            },
            timeout=30,
        )
        resp.raise_for_status()
        return resp.json()


class LocalReductionBackend:
    def submit(self, args: List[str]) -> None:
        parsed = ql.QL.get_parser().parse_args(args)
        ql.QL.main(parsed)


class RemoteReductionBackend:
    """Reduction backend that submits jobs to a remote :mod:`HTTPserver`."""

    def __init__(self, base_url: str) -> None:
        self.base_url = base_url.rstrip("/")

    def submit(self, args: List[str]) -> None:
        """Fire a reduction job on the remote server.

        The server runs the job on a daemon thread and returns immediately with
        a job ID.  Since the plugin already wraps ``submit`` in its own daemon
        thread and polls for output files via the filesystem, we do not need to
        track the job ID here — the fire-and-forget pattern is sufficient.
        """
        import requests
        resp = requests.post(
            f"{self.base_url}/api/submit",
            json={"args": args},
            timeout=30,
        )
        resp.raise_for_status()
