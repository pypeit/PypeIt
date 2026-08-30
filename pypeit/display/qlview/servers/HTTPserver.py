"""
PypeIt QuickLook HTTP server.

Exposes the file-browser and reduction backends over a simple REST API so
that a remote Ginga instance (or any HTTP client) can browse files and
launch PypeIt quicklook reductions on the machine where raw data and
calibrations reside.

Endpoints
---------
GET  /api/health
    Liveness probe — exempted from API-key auth.
    Returns ``{"status": "ok"}``.

GET  /api/list_dir?path=<path>
    List a directory.  Returns ``{"files": ["/abs/path1", ...]}``.

GET  /api/stat_info?path=<path>
    Stat a path.  Returns serialised stat fields.

GET  /api/header_info?path=<path>&instrument=<name>&mode=<raw|reduced>
    Read FITS or calibration-directory metadata.

GET  /api/check_log?path=<log_path>&failure_string=<str>
    Returns ``{"failed": true|false}``.

GET  /api/glob?dir=<dir>&pattern=<pattern>
    Glob a server-side directory.  Returns ``{"files": [...]}``.

POST /api/submit   body: {"args": ["keck_deimos", "--raw_files", ...]}
    Submit a reduction.  Returns ``{"job_id": "<uuid>", "status": "running"}``.

GET  /api/job/<job_id>
    Poll job status.  Returns status/error/completed_at fields.

DELETE /api/job/<job_id>
    Remove a completed or failed job from the store.

Security
--------
- Bind to ``127.0.0.1`` by default (local-only).  Pass ``--host 0.0.0.0``
  to listen on all interfaces.
- Use ``--api-key <token>`` to require ``Authorization: Bearer <token>`` on
  every request (except ``/api/health``).
- Use ``--allowed-roots /data/raw /data/cal`` to restrict path arguments to
  specific filesystem subtrees.  Without this flag all local paths are
  accessible.

Production use
--------------
For concurrent workloads, run behind a WSGI server::

    gunicorn -w 4 pypeit.display.qlview.servers.HTTPserver:app

Usage
-----
::

    python -m pypeit.display.qlview.servers.HTTPserver \\
        --host 127.0.0.1 --port 5000 \\
        --api-key secret123 \\
        --allowed-roots /data/raw /data/cal /data/redux
"""

from __future__ import annotations

import argparse
import logging
import os
import threading
import time
import uuid
from typing import Any, Dict, List

from flask import Flask, jsonify, request

from pypeit.display.qlview.backends import LocalFileBrowserBackend, LocalReductionBackend
from pypeit.display.qlview.instruments import InstrumentRegistry

# ---------------------------------------------------------------------------
# Flask application
# ---------------------------------------------------------------------------

app = Flask(__name__)

# ---------------------------------------------------------------------------
# Module-level singletons (configured at startup via run())
# ---------------------------------------------------------------------------

_file_backend = LocalFileBrowserBackend()
_reduction_backend = LocalReductionBackend()
_registry: InstrumentRegistry | None = None
_registry_lock = threading.Lock()

_api_key: str | None = None
_allowed_roots: List[str] = []  # empty = no restriction

# Job store: job_id -> {"status", "error", "completed_at"}
_jobs: Dict[str, Dict[str, Any]] = {}
_jobs_lock = threading.Lock()

_JOB_TTL_SECONDS = 3600  # evict completed/failed jobs after 1 hour


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_registry() -> InstrumentRegistry:
    global _registry
    with _registry_lock:
        if _registry is None:
            _registry = InstrumentRegistry(app.logger)
    return _registry


def _validate_path(path: str) -> str:
    """Resolve *path* and verify it is under an allowed root.

    Returns the resolved absolute path.  Raises :class:`PermissionError` if
    ``_allowed_roots`` is non-empty and the path is not under any of them.
    """
    resolved = os.path.realpath(path)
    if _allowed_roots:
        if not any(resolved.startswith(root) for root in _allowed_roots):
            raise PermissionError(f"Path not allowed: {path!r}")
    return resolved


@app.before_request
def _require_auth():
    """Reject unauthenticated requests when an API key is configured.

    ``/api/health`` is exempt so that liveness probes work without credentials.
    """
    if _api_key is None:
        return  # auth disabled
    if request.endpoint == "health":
        return  # exempt
    auth = request.headers.get("Authorization", "")
    if auth != f"Bearer {_api_key}":
        return jsonify({"error": "Unauthorized"}), 401


def _evict_old_jobs() -> None:
    """Background thread: remove completed/failed jobs older than the TTL."""
    while True:
        time.sleep(60)
        cutoff = time.time() - _JOB_TTL_SECONDS
        with _jobs_lock:
            stale = [
                jid for jid, job in _jobs.items()
                if job["completed_at"] is not None and job["completed_at"] < cutoff
            ]
            for jid in stale:
                del _jobs[jid]


# Start the eviction thread once at module load (daemon so it doesn't block exit).
threading.Thread(target=_evict_old_jobs, daemon=True, name="ql-job-evict").start()


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

@app.route("/api/health")
def health():
    """Liveness probe (no auth required)."""
    return jsonify({"status": "ok"})


@app.route("/api/list_dir")
def list_dir():
    path = request.args.get("path", "")
    if not path:
        return jsonify({"error": "Missing required query parameter: path"}), 400
    try:
        safe = _validate_path(path)
        files = _file_backend.list_dir(safe)
        return jsonify({"files": files})
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/stat_info")
def stat_info():
    path = request.args.get("path", "")
    if not path:
        return jsonify({"error": "Missing required query parameter: path"}), 400
    try:
        safe = _validate_path(path)
        s = _file_backend.stat_info(safe)
        return jsonify({
            "st_mode": s.st_mode,
            "st_size": s.st_size,
            "st_mtime": s.st_mtime,
            "st_uid": s.st_uid,
            "st_gid": s.st_gid,
        })
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except OSError as exc:
        return jsonify({"error": str(exc)}), 404
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/header_info")
def header_info():
    path = request.args.get("path", "")
    instrument_name = request.args.get("instrument", "DEIMOS")
    mode = request.args.get("mode", "raw")

    if not path:
        return jsonify({"error": "Missing required query parameter: path"}), 400

    try:
        safe = _validate_path(path)
        instrument = _get_registry().create(instrument_name)
        info = _file_backend.get_header_info(safe, instrument, mode=mode)
        serialisable = {
            k: (str(v) if v is not None else "N/A")
            for k, v in info.items()
        }
        return jsonify(serialisable)
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/check_log")
def check_log():
    """Check a log file for a failure string.

    Query parameters
    ----------------
    path : str
        Path to the log file on the server filesystem.
    failure_string : str
        Substring to search for.
    """
    path = request.args.get("path", "")
    failure_string = request.args.get("failure_string", "")
    if not path:
        return jsonify({"error": "Missing required query parameter: path"}), 400
    if not failure_string:
        return jsonify({"error": "Missing required query parameter: failure_string"}), 400
    try:
        safe = _validate_path(path)
        failed = _file_backend.check_log_for_failure(safe, failure_string)
        return jsonify({"failed": failed})
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/glob")
def glob_dir():
    """Glob a server-side directory.

    Query parameters
    ----------------
    dir : str
        Directory to search.
    pattern : str
        Glob pattern (e.g. ``spec1d*.fits*``).
    """
    dir_path = request.args.get("dir", "")
    pattern = request.args.get("pattern", "")
    if not dir_path:
        return jsonify({"error": "Missing required query parameter: dir"}), 400
    if not pattern:
        return jsonify({"error": "Missing required query parameter: pattern"}), 400
    try:
        safe = _validate_path(dir_path)
        files = _file_backend.glob(safe, pattern)
        return jsonify({"files": files})
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/submit", methods=["POST"])
def submit():
    """Submit a PypeIt QuickLook reduction.

    Request body (JSON)
    -------------------
    args : list[str]
        Command-line arguments passed to ``ql.QL``.
    """
    data = request.get_json(force=True, silent=True)
    if not data or "args" not in data:
        return jsonify({"error": "Request body must be JSON with an 'args' key"}), 400

    args = data["args"]
    if not isinstance(args, list) or not all(isinstance(a, str) for a in args):
        return jsonify({"error": "'args' must be a list of strings"}), 400

    # Validate any path-like arguments against allowed roots.
    path_flags = {"--redux_path", "--raw_path", "--setup_calib_dir", "--log_file"}
    for flag in path_flags:
        if flag in args:
            idx = args.index(flag)
            if idx + 1 < len(args):
                try:
                    _validate_path(args[idx + 1])
                except PermissionError as exc:
                    return jsonify({"error": str(exc)}), 403

    job_id = str(uuid.uuid4())
    with _jobs_lock:
        _jobs[job_id] = {"status": "running", "error": None, "completed_at": None}

    def _run() -> None:
        try:
            _reduction_backend.submit(args)
            with _jobs_lock:
                _jobs[job_id]["status"] = "completed"
                _jobs[job_id]["completed_at"] = time.time()
        except Exception as exc:
            with _jobs_lock:
                _jobs[job_id]["status"] = "failed"
                _jobs[job_id]["error"] = str(exc)
                _jobs[job_id]["completed_at"] = time.time()

    threading.Thread(target=_run, daemon=True, name=f"ql-job-{job_id[:8]}").start()
    return jsonify({"job_id": job_id, "status": "running"}), 202


@app.route("/api/job/<job_id>")
def job_status(job_id: str):
    """Poll the status of a submitted reduction job."""
    with _jobs_lock:
        job = _jobs.get(job_id)
    if job is None:
        return jsonify({"error": f"Job not found: {job_id}"}), 404
    return jsonify({"job_id": job_id, **job})


@app.route("/api/job/<job_id>", methods=["DELETE"])
def delete_job(job_id: str):
    """Remove a job from the store."""
    with _jobs_lock:
        if job_id not in _jobs:
            return jsonify({"error": f"Job not found: {job_id}"}), 404
        del _jobs[job_id]
    return "", 204


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run(
    host: str = "127.0.0.1",
    port: int = 5000,
    debug: bool = False,
    threaded: bool = True,
    api_key: str | None = None,
    allowed_roots: List[str] | None = None,
) -> None:
    """Start the Flask development server.

    Parameters
    ----------
    host : str
        Interface to bind.  Defaults to ``127.0.0.1`` (local-only).
        Pass ``"0.0.0.0"`` to listen on all interfaces.
    port : int
        TCP port number.
    debug : bool
        Enable Flask debug mode (auto-reload, verbose errors).
        Do **not** use in production.
    threaded : bool
        Enable Werkzeug threaded mode so concurrent requests are handled
        without queuing.  Defaults to ``True``.
        For production use, run under gunicorn instead::

            gunicorn -w 4 pypeit.display.qlview.servers.HTTPserver:app
    api_key : str or None
        Shared secret token.  When set, every request (except ``/api/health``)
        must include ``Authorization: Bearer <token>``.
    allowed_roots : list of str or None
        Filesystem roots that path arguments are restricted to.  An empty list
        (the default) disables path sandboxing.
    """
    global _api_key, _allowed_roots

    app.logger.setLevel(logging.INFO)

    _api_key = api_key or None
    _allowed_roots = [os.path.realpath(r) for r in (allowed_roots or [])]

    app.run(host=host, port=port, debug=debug, threaded=threaded)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="PypeIt QuickLook HTTP backend server"
    )
    parser.add_argument(
        "--host", default="127.0.0.1",
        help="Bind address (default: 127.0.0.1 — local only)",
    )
    parser.add_argument("--port", type=int, default=5000, help="Port number (default: 5000)")
    parser.add_argument("--debug", action="store_true", help="Enable Flask debug mode")
    parser.add_argument(
        "--no-threaded", dest="threaded", action="store_false", default=True,
        help="Disable Werkzeug threaded mode (useful for debugging)",
    )
    parser.add_argument(
        "--api-key", default=None, dest="api_key",
        help="Shared secret token for Bearer auth (no auth if omitted)",
    )
    parser.add_argument(
        "--allowed-roots", nargs="*", default=[], dest="allowed_roots",
        metavar="PATH",
        help="Filesystem roots that path arguments must be under (no restriction if omitted)",
    )
    parsed = parser.parse_args()
    run(
        host=parsed.host,
        port=parsed.port,
        debug=parsed.debug,
        threaded=parsed.threaded,
        api_key=parsed.api_key,
        allowed_roots=parsed.allowed_roots,
    )
