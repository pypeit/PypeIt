from __future__ import annotations

import glob
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Protocol, TYPE_CHECKING

if TYPE_CHECKING:
    from .instruments import Instrument

from pypeit.scripts import ql

# requests is optional — only needed for the remote backends.
try:
    import requests as _requests
except ImportError:
    _requests = None


# ---------------------------------------------------------------------------
# Shared types
# ---------------------------------------------------------------------------

@dataclass
class StatResult:
    """Filesystem stat fields consumed by :class:`FileBrowserController`."""
    st_mode: int
    st_size: int
    st_mtime: float
    st_uid: int
    st_gid: int


# ---------------------------------------------------------------------------
# Backend protocols
# ---------------------------------------------------------------------------

class FileBrowserBackend(Protocol):
    def list_dir(self, path: str) -> List[str]:
        ...

    def stat_info(self, path: str) -> StatResult:
        ...

    def get_header_info(self, path: str, instrument: "Instrument", mode: str = "raw") -> Dict[str, object]:
        ...

    def check_log_for_failure(self, log_path: str, failure_string: str) -> bool:
        ...

    def glob(self, dir_path: str, pattern: str) -> List[str]:
        ...


class ReductionBackend(Protocol):
    def submit(self, args: List[str]) -> None:
        ...


# ---------------------------------------------------------------------------
# Local implementations
# ---------------------------------------------------------------------------

class LocalFileBrowserBackend:
    def list_dir(self, path: str) -> List[str]:
        fullpath = os.path.join(os.path.abspath(path), "*")
        filelist = list(glob.glob(fullpath))
        filelist.sort(key=lambda s: s.lower())
        filelist.insert(0, os.path.join(os.path.abspath(path), ".."))
        return filelist

    def stat_info(self, path: str) -> StatResult:
        s = os.stat(path)
        return StatResult(
            st_mode=s.st_mode,
            st_size=s.st_size,
            st_mtime=s.st_mtime,
            st_uid=s.st_uid,
            st_gid=s.st_gid,
        )

    def get_header_info(self, path: str, instrument: "Instrument", mode: str = "raw") -> Dict[str, object]:
        if mode == "reduced":
            return instrument.get_reduced_info(path)
        return instrument.get_raw_info(path)

    def check_log_for_failure(self, log_path: str, failure_string: str) -> bool:
        try:
            with open(log_path) as fh:
                return failure_string in fh.read()
        except OSError:
            return False

    def glob(self, dir_path: str, pattern: str) -> List[str]:
        p = Path(dir_path)
        if not p.exists():
            return []
        return sorted(str(f) for f in p.glob(pattern))


# ---------------------------------------------------------------------------
# Remote implementations
# ---------------------------------------------------------------------------

def _raise_for_status(resp) -> None:
    """Raise a :class:`RuntimeError` with the server's JSON error body, if any."""
    try:
        resp.raise_for_status()
    except _requests.HTTPError as exc:
        try:
            msg = resp.json().get("error", str(exc))
        except ValueError:
            msg = str(exc)
        raise RuntimeError(f"Server error ({resp.status_code}): {msg}") from exc


class RemoteFileBrowserBackend:
    """File-browser backend that delegates to a remote :mod:`HTTPserver`."""

    def __init__(self, base_url: str, api_key: str | None = None) -> None:
        if _requests is None:
            raise RuntimeError("The 'requests' package is required for remote backends.")
        self.base_url = base_url.rstrip("/")
        self.api_key = api_key

    def _headers(self) -> Dict[str, str]:
        if self.api_key:
            return {"Authorization": f"Bearer {self.api_key}"}
        return {}

    def list_dir(self, path: str) -> List[str]:
        resp = _requests.get(
            f"{self.base_url}/api/list_dir",
            params={"path": path},
            headers=self._headers(),
            timeout=30,
        )
        _raise_for_status(resp)
        return resp.json()["files"]

    def stat_info(self, path: str) -> StatResult:
        resp = _requests.get(
            f"{self.base_url}/api/stat_info",
            params={"path": path},
            headers=self._headers(),
            timeout=30,
        )
        _raise_for_status(resp)
        data = resp.json()
        return StatResult(
            st_mode=data["st_mode"],
            st_size=data["st_size"],
            st_mtime=data["st_mtime"],
            st_uid=data.get("st_uid", 0),
            st_gid=data.get("st_gid", 0),
        )

    def get_header_info(self, path: str, instrument: "Instrument", mode: str = "raw") -> Dict[str, object]:
        resp = _requests.get(
            f"{self.base_url}/api/header_info",
            params={
                "path": path,
                "instrument": instrument.__class__.__name__,
                "mode": mode,
            },
            headers=self._headers(),
            timeout=30,
        )
        _raise_for_status(resp)
        return resp.json()

    def check_log_for_failure(self, log_path: str, failure_string: str) -> bool:
        resp = _requests.get(
            f"{self.base_url}/api/check_log",
            params={"path": log_path, "failure_string": failure_string},
            headers=self._headers(),
            timeout=30,
        )
        _raise_for_status(resp)
        return resp.json()["failed"]

    def glob(self, dir_path: str, pattern: str) -> List[str]:
        resp = _requests.get(
            f"{self.base_url}/api/glob",
            params={"dir": dir_path, "pattern": pattern},
            headers=self._headers(),
            timeout=30,
        )
        _raise_for_status(resp)
        return resp.json()["files"]


class LocalReductionBackend:
    def submit(self, args: List[str]) -> None:
        parsed = ql.QL.get_parser().parse_args(args)
        ql.QL.main(parsed)


class RemoteReductionBackend:
    """Reduction backend that submits jobs to a remote :mod:`HTTPserver`."""

    def __init__(self, base_url: str, api_key: str | None = None) -> None:
        if _requests is None:
            raise RuntimeError("The 'requests' package is required for remote backends.")
        self.base_url = base_url.rstrip("/")
        self.api_key = api_key

    def _headers(self) -> Dict[str, str]:
        if self.api_key:
            return {"Authorization": f"Bearer {self.api_key}"}
        return {}

    def submit(self, args: List[str]) -> None:
        """Fire a reduction job on the remote server.

        The server runs the job on a daemon thread and returns immediately with
        a job ID.  Since the plugin already wraps ``submit`` in its own daemon
        thread and polls for output files via the filesystem backend, we do not
        need to track the job ID here.
        """
        resp = _requests.post(
            f"{self.base_url}/api/submit",
            json={"args": args},
            headers=self._headers(),
            timeout=30,
        )
        _raise_for_status(resp)
