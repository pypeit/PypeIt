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
        """Return a sorted listing of all entries under *path*.

        Parameters
        ----------
        path : str
            Directory path to list.

        Returns
        -------
        list of str
            Absolute paths of directory entries, with a parent-directory entry
            (``".."``) prepended.
        """
        ...

    def stat_info(self, path: str) -> StatResult:
        """Return filesystem stat fields for *path*.

        Parameters
        ----------
        path : str
            Path to stat.

        Returns
        -------
        StatResult
            Stat fields (mode, size, mtime, uid, gid) for the given path.
        """
        ...

    def get_header_info(self, path: str, instrument: "Instrument", mode: str = "raw") -> Dict[str, object]:
        """Return header metadata for a FITS file.

        Parameters
        ----------
        path : str
            Path to the FITS file.
        instrument : Instrument
            Instrument object used to extract metadata.
        mode : str, optional
            ``"raw"`` (default) to read raw-frame headers; ``"reduced"`` to
            read reduced-product headers.

        Returns
        -------
        dict
            Mapping of header keyword names to their values.
        """
        ...

    def check_log_for_failure(self, log_path: str, failure_string: str) -> bool:
        """Check whether a log file contains a failure indicator string.

        Parameters
        ----------
        log_path : str
            Path to the log file.
        failure_string : str
            String whose presence indicates a failure.

        Returns
        -------
        bool
            ``True`` if *failure_string* is found in the log; ``False``
            otherwise (including when the file cannot be read).
        """
        ...

    def glob(self, dir_path: str, pattern: str) -> List[str]:
        """Return sorted paths matching *pattern* inside *dir_path*.

        Parameters
        ----------
        dir_path : str
            Directory to search.
        pattern : str
            Glob pattern relative to *dir_path*.

        Returns
        -------
        list of str
            Sorted absolute paths of all matches.  Returns an empty list when
            *dir_path* does not exist.
        """
        ...


class ReductionBackend(Protocol):
    def submit(self, args: List[str]) -> None:
        """Submit a quicklook reduction job.

        Parameters
        ----------
        args : list of str
            Command-line argument list passed to the quicklook reduction
            script (``pypeit_ql``).
        """
        ...


# ---------------------------------------------------------------------------
# Local implementations
# ---------------------------------------------------------------------------

class LocalFileBrowserBackend:
    def list_dir(self, path: str) -> List[str]:
        """Implements :meth:`FileBrowserBackend.list_dir` for the local filesystem.

        Parameters
        ----------
        path : str
            Directory path to list.

        Returns
        -------
        list of str
            Sorted absolute paths of directory entries with a ``".."`` parent
            entry prepended.
        """
        fullpath = os.path.join(os.path.abspath(path), "*")
        filelist = list(glob.glob(fullpath))
        filelist.sort(key=lambda s: s.lower())
        filelist.insert(0, os.path.join(os.path.abspath(path), ".."))
        return filelist

    def stat_info(self, path: str) -> StatResult:
        """Implements :meth:`FileBrowserBackend.stat_info` for the local filesystem.

        Parameters
        ----------
        path : str
            Path to stat.

        Returns
        -------
        StatResult
            Stat fields obtained via :func:`os.stat`.
        """
        s = os.stat(path)
        return StatResult(
            st_mode=s.st_mode,
            st_size=s.st_size,
            st_mtime=s.st_mtime,
            st_uid=s.st_uid,
            st_gid=s.st_gid,
        )

    def get_header_info(self, path: str, instrument: "Instrument", mode: str = "raw") -> Dict[str, object]:
        """Implements :meth:`FileBrowserBackend.get_header_info` for the local filesystem.

        Parameters
        ----------
        path : str
            Path to the FITS file.
        instrument : Instrument
            Instrument object used to extract metadata.
        mode : str, optional
            ``"raw"`` (default) delegates to
            :meth:`~pypeit.display.qlview.instruments.Instrument.get_raw_info`;
            ``"reduced"`` delegates to
            :meth:`~pypeit.display.qlview.instruments.Instrument.get_reduced_info`.

        Returns
        -------
        dict
            Mapping of header keyword names to their values.
        """
        if mode == "reduced":
            return instrument.get_reduced_info(path)
        return instrument.get_raw_info(path)

    def check_log_for_failure(self, log_path: str, failure_string: str) -> bool:
        """Implements :meth:`FileBrowserBackend.check_log_for_failure` for the local filesystem.

        Parameters
        ----------
        log_path : str
            Path to the log file.
        failure_string : str
            String whose presence indicates a failure.

        Returns
        -------
        bool
            ``True`` if *failure_string* is found; ``False`` if the string is
            absent or the file cannot be opened.
        """
        try:
            with open(log_path) as fh:
                return failure_string in fh.read()
        except OSError:
            return False

    def glob(self, dir_path: str, pattern: str) -> List[str]:
        """Implements :meth:`FileBrowserBackend.glob` for the local filesystem.

        Parameters
        ----------
        dir_path : str
            Directory to search.
        pattern : str
            Glob pattern relative to *dir_path*.

        Returns
        -------
        list of str
            Sorted absolute paths of all matches.  Returns an empty list when
            *dir_path* does not exist.
        """
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
        """
        Parameters
        ----------
        base_url : str
            Base URL of the remote HTTP server (e.g. ``"http://host:8765"``).
            Trailing slashes are stripped automatically.
        api_key : str, optional
            API key sent as a Bearer token in the ``Authorization`` header.
            Pass ``None`` (default) for unauthenticated servers.

        Raises
        ------
        RuntimeError
            If the ``requests`` package is not installed.
        """
        if _requests is None:
            raise RuntimeError("The 'requests' package is required for remote backends.")
        self.base_url = base_url.rstrip("/")
        self.api_key = api_key

    def _headers(self) -> Dict[str, str]:
        """Return the Authorization header dict used for all HTTP requests.

        Returns
        -------
        dict
            ``{"Authorization": "Bearer <api_key>"}`` when an API key is set;
            empty dict otherwise.
        """
        if self.api_key:
            return {"Authorization": f"Bearer {self.api_key}"}
        return {}

    def list_dir(self, path: str) -> List[str]:
        """Implements :meth:`FileBrowserBackend.list_dir` over HTTP.

        Sends a GET request to ``<base_url>/api/list_dir``.

        Parameters
        ----------
        path : str
            Directory path on the remote host.

        Returns
        -------
        list of str
            Directory entries returned by the server.

        Raises
        ------
        RuntimeError
            If the server returns a non-2xx HTTP status.
        """
        resp = _requests.get(
            f"{self.base_url}/api/list_dir",
            params={"path": path},
            headers=self._headers(),
            timeout=30,
        )
        _raise_for_status(resp)
        return resp.json()["files"]

    def stat_info(self, path: str) -> StatResult:
        """Implements :meth:`FileBrowserBackend.stat_info` over HTTP.

        Sends a GET request to ``<base_url>/api/stat_info``.

        Parameters
        ----------
        path : str
            Path to stat on the remote host.

        Returns
        -------
        StatResult
            Stat fields returned by the server.

        Raises
        ------
        RuntimeError
            If the server returns a non-2xx HTTP status.
        """
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
        """Implements :meth:`FileBrowserBackend.get_header_info` over HTTP.

        Sends a GET request to ``<base_url>/api/header_info``.  The instrument
        class name is passed as a query parameter so the server can instantiate
        the correct instrument.

        Parameters
        ----------
        path : str
            Path to the FITS file on the remote host.
        instrument : Instrument
            Instrument whose class name is forwarded to the server.
        mode : str, optional
            ``"raw"`` (default) or ``"reduced"``.

        Returns
        -------
        dict
            Mapping of header keyword names to their values.

        Raises
        ------
        RuntimeError
            If the server returns a non-2xx HTTP status.
        """
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
        """Implements :meth:`FileBrowserBackend.check_log_for_failure` over HTTP.

        Sends a GET request to ``<base_url>/api/check_log``.

        Parameters
        ----------
        log_path : str
            Path to the log file on the remote host.
        failure_string : str
            String whose presence indicates a failure.

        Returns
        -------
        bool
            ``True`` if *failure_string* is found according to the server.

        Raises
        ------
        RuntimeError
            If the server returns a non-2xx HTTP status.
        """
        resp = _requests.get(
            f"{self.base_url}/api/check_log",
            params={"path": log_path, "failure_string": failure_string},
            headers=self._headers(),
            timeout=30,
        )
        _raise_for_status(resp)
        return resp.json()["failed"]

    def glob(self, dir_path: str, pattern: str) -> List[str]:
        """Implements :meth:`FileBrowserBackend.glob` over HTTP.

        Sends a GET request to ``<base_url>/api/glob``.

        Parameters
        ----------
        dir_path : str
            Directory to search on the remote host.
        pattern : str
            Glob pattern relative to *dir_path*.

        Returns
        -------
        list of str
            Matching paths returned by the server.

        Raises
        ------
        RuntimeError
            If the server returns a non-2xx HTTP status.
        """
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
        """Implements :meth:`ReductionBackend.submit` synchronously on the calling thread.

        Parses *args* with the quicklook argument parser and calls
        :meth:`~pypeit.scripts.ql.QL.main` directly.  Because this runs
        synchronously on the calling thread, callers that need non-blocking
        behaviour should invoke this method from a background thread.

        Parameters
        ----------
        args : list of str
            Command-line argument list accepted by the ``pypeit_ql`` script.
        """
        parsed = ql.QL.get_parser().parse_args(args)
        ql.QL.main(parsed)


class RemoteReductionBackend:
    """Reduction backend that submits jobs to a remote :mod:`HTTPserver`."""

    def __init__(self, base_url: str, api_key: str | None = None) -> None:
        """
        Parameters
        ----------
        base_url : str
            Base URL of the remote HTTP server (e.g. ``"http://host:8765"``).
            Trailing slashes are stripped automatically.
        api_key : str, optional
            API key sent as a Bearer token in the ``Authorization`` header.
            Pass ``None`` (default) for unauthenticated servers.

        Raises
        ------
        RuntimeError
            If the ``requests`` package is not installed.
        """
        if _requests is None:
            raise RuntimeError("The 'requests' package is required for remote backends.")
        self.base_url = base_url.rstrip("/")
        self.api_key = api_key

    def _headers(self) -> Dict[str, str]:
        """Return the Authorization header dict used for all HTTP requests.

        Returns
        -------
        dict
            ``{"Authorization": "Bearer <api_key>"}`` when an API key is set;
            empty dict otherwise.
        """
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
