"""
PypeIt QuickLook HTTP server.

Exposes the file-browser and reduction backends over a simple REST API so
that a remote Ginga instance (or any HTTP client) can browse files and
launch PypeIt quicklook reductions on the machine where raw data and
calibrations reside.

Endpoints
---------
GET  /api/health
    Simple liveness probe.  Returns ``{"status": "ok"}``.

GET  /api/list_dir?path=<path>
    List a directory.  Returns ``{"files": ["/abs/path1", ...]}``.

GET  /api/stat_info?path=<path>
    Stat a path.  Returns serialised ``os.stat_result`` fields.

GET  /api/header_info?path=<path>&instrument=<name>&mode=<raw|reduced>
    Read FITS or calibration-directory metadata.
    *instrument* must match a key in ``InstrumentRegistry`` (e.g. ``DEIMOS``).

POST /api/submit   body: {"args": ["keck_deimos", "--raw_files", ...]}
    Submit a reduction.  The job runs on a daemon thread.
    Returns ``{"job_id": "<uuid>", "status": "running"}``.

GET  /api/job/<job_id>
    Poll a running or finished job.
    Returns ``{"job_id": ..., "status": "running"|"completed"|"failed",
                "error": null|"<message>"}``.

Usage
-----
Run directly::

    python -m pypeit.display.qlview.servers.HTTPserver --port 5000

Or import and call :func:`run`::

    from pypeit.display.qlview.servers.HTTPserver import run
    run(host="0.0.0.0", port=5000)
"""

from __future__ import annotations

import argparse
import logging
import threading
import uuid
from typing import Any, Dict

from flask import Flask, jsonify, request

from pypeit.display.qlview.backends import LocalFileBrowserBackend, LocalReductionBackend
from pypeit.display.qlview.instruments import InstrumentRegistry

# ---------------------------------------------------------------------------
# Flask application
# ---------------------------------------------------------------------------

app = Flask(__name__)

# ---------------------------------------------------------------------------
# Module-level singletons
# ---------------------------------------------------------------------------

_file_backend = LocalFileBrowserBackend()
_reduction_backend = LocalReductionBackend()
_registry: InstrumentRegistry | None = None
_registry_lock = threading.Lock()

# Job store: job_id -> {"status": "running"|"completed"|"failed", "error": str|None}
_jobs: Dict[str, Dict[str, Any]] = {}
_jobs_lock = threading.Lock()


def _get_registry() -> InstrumentRegistry:
    global _registry
    with _registry_lock:
        if _registry is None:
            _registry = InstrumentRegistry(logging.getLogger("pypeit"))
    return _registry


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

@app.route("/api/health")
def health():
    """Liveness probe."""
    return jsonify({"status": "ok"})


@app.route("/api/list_dir")
def list_dir():
    """
    List a directory.

    Query parameters
    ----------------
    path : str
        Absolute or relative path to the directory to list.
    """
    path = request.args.get("path", "")
    if not path:
        return jsonify({"error": "Missing required query parameter: path"}), 400
    try:
        files = _file_backend.list_dir(path)
        return jsonify({"files": files})
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/stat_info")
def stat_info():
    """
    Stat a filesystem path.

    Query parameters
    ----------------
    path : str
        Path to stat.
    """
    path = request.args.get("path", "")
    if not path:
        return jsonify({"error": "Missing required query parameter: path"}), 400
    try:
        s = _file_backend.stat_info(path)
        return jsonify({
            "st_mode": s.st_mode,
            "st_size": s.st_size,
            "st_mtime": s.st_mtime,
            "st_uid": s.st_uid,
            "st_gid": s.st_gid,
        })
    except OSError as exc:
        return jsonify({"error": str(exc)}), 404
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/header_info")
def header_info():
    """
    Read FITS or calibration-directory metadata.

    Query parameters
    ----------------
    path : str
        Path to the file or directory.
    instrument : str
        Instrument registry key (e.g. ``DEIMOS``, ``MOSFIRE``).
    mode : str
        ``"raw"`` (default) or ``"reduced"``.
    """
    path = request.args.get("path", "")
    instrument_name = request.args.get("instrument", "DEIMOS")
    mode = request.args.get("mode", "raw")

    if not path:
        return jsonify({"error": "Missing required query parameter: path"}), 400

    try:
        instrument = _get_registry().create(instrument_name)
        info = _file_backend.get_header_info(path, instrument, mode=mode)
        # Ensure every value is JSON-serialisable
        serialisable = {
            k: (str(v) if v is not None else "N/A")
            for k, v in info.items()
        }
        return jsonify(serialisable)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/submit", methods=["POST"])
def submit():
    """
    Submit a PypeIt QuickLook reduction.

    Request body (JSON)
    -------------------
    args : list[str]
        Command-line arguments to pass to ``ql.QL``, e.g.
        ``["keck_deimos", "--raw_files", "d1234.fits", ...]``.

    Returns
    -------
    JSON with ``job_id`` and initial ``status`` of ``"running"``.
    """
    data = request.get_json(force=True, silent=True)
    if not data or "args" not in data:
        return jsonify({"error": "Request body must be JSON with an 'args' key"}), 400

    args = data["args"]
    if not isinstance(args, list) or not all(isinstance(a, str) for a in args):
        return jsonify({"error": "'args' must be a list of strings"}), 400

    job_id = str(uuid.uuid4())
    with _jobs_lock:
        _jobs[job_id] = {"status": "running", "error": None}

    def _run() -> None:
        try:
            _reduction_backend.submit(args)
            with _jobs_lock:
                _jobs[job_id]["status"] = "completed"
        except Exception as exc:
            with _jobs_lock:
                _jobs[job_id]["status"] = "failed"
                _jobs[job_id]["error"] = str(exc)

    threading.Thread(target=_run, daemon=True, name=f"ql-job-{job_id[:8]}").start()
    return jsonify({"job_id": job_id, "status": "running"}), 202


@app.route("/api/job/<job_id>")
def job_status(job_id: str):
    """
    Poll the status of a submitted reduction job.

    Path parameters
    ---------------
    job_id : str
        UUID returned by ``POST /api/submit``.
    """
    with _jobs_lock:
        job = _jobs.get(job_id)
    if job is None:
        return jsonify({"error": f"Job not found: {job_id}"}), 404
    return jsonify({"job_id": job_id, **job})


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run(host: str = "0.0.0.0", port: int = 5000, debug: bool = False) -> None:
    """Start the Flask development server.

    Parameters
    ----------
    host : str
        Interface to bind (``"0.0.0.0"`` listens on all interfaces).
    port : int
        TCP port number.
    debug : bool
        Enable Flask debug mode (auto-reload, verbose errors).  Do **not**
        use in production.
    """
    app.run(host=host, port=port, debug=debug)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="PypeIt QuickLook HTTP backend server"
    )
    parser.add_argument("--host", default="0.0.0.0", help="Bind address (default: 0.0.0.0)")
    parser.add_argument("--port", type=int, default=5000, help="Port number (default: 5000)")
    parser.add_argument("--debug", action="store_true", help="Enable Flask debug mode")
    parsed = parser.parse_args()
    run(host=parsed.host, port=parsed.port, debug=parsed.debug)
