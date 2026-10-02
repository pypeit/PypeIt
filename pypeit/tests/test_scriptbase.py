"""Tests for shared command-line script infrastructure."""

import matplotlib

from pypeit.scripts import scriptbase


def test_configure_matplotlib_headless(monkeypatch):
    calls = []
    monkeypatch.setattr(matplotlib, 'use', lambda backend, force=False: calls.append((backend, force)))

    scriptbase.configure_matplotlib(show=False)

    assert calls == [('Agg', True)]


def test_configure_matplotlib_show_preserves_backend(monkeypatch):
    calls = []
    monkeypatch.setattr(matplotlib, 'use', lambda backend, force=False: calls.append((backend, force)))

    scriptbase.configure_matplotlib(show=True)

    assert calls == []
