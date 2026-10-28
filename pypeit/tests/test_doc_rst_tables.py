"""Check the AstroPy-formatted tables produced by the documentation scripts."""

import importlib
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

from docutils import nodes
from docutils.core import publish_doctree


DOC_SCRIPTS = Path(__file__).resolve().parents[2] / 'doc' / 'scripts'
PROJECT_FILE = DOC_SCRIPTS.parents[1] / 'pyproject.toml'


def load_script(monkeypatch, name):
    monkeypatch.syspath_prepend(str(DOC_SCRIPTS))
    return importlib.import_module(name)


def table_rows(source):
    """Return the rendered text of every cell in each RST table."""
    warnings = StringIO()
    tree = publish_doctree(source, settings_overrides={
        'warning_stream': warnings, 'report_level': 2, 'halt_level': 6})
    assert 'Malformed table' not in warnings.getvalue()
    return [
        [[entry.astext() for entry in row.findall(nodes.entry)]
         for row in table.findall(nodes.row)]
        for table in tree.findall(nodes.table)
    ]


def test_atmospheric_extinction_table(monkeypatch, tmp_path):
    script = load_script(monkeypatch, 'build_atmexttbl_rst')
    script.write_tables(tmp_path)
    rows = table_rows((tmp_path / 'atmext_table.rst').read_text())
    assert len(rows) == 1
    assert len(rows[0]) > 1
    assert 'File' in rows[0][0]


def test_bitmask_table(monkeypatch, tmp_path):
    script = load_script(monkeypatch, 'build_bitmask_rst')
    from pypeit.images.imagebitmask import ImageBitMask

    script.write_bitmask_table(ImageBitMask, tmp_path)
    rows = table_rows((tmp_path / 'imagebitmask_table.rst').read_text())
    assert len(rows) == 1
    assert rows[0][0] == ['Bit Name', 'Bit Number', 'Decimal Value', 'Description']
    assert {row[0] for row in rows[0][1:]} == set(ImageBitMask().bits)


def test_cache_data_table(monkeypatch, tmp_path):
    script = load_script(monkeypatch, 'build_cache_data_tbl')
    from pypeit.pypeitdata import PypeItDataPaths

    output = tmp_path / 'data_dir.rst'
    script.write_table(output)
    rows = table_rows(output.read_text())
    assert len(rows) == 1
    assert rows[0][0] == ['Reference', 'Subdirectory', 'Host']
    assert [row[0] for row in rows[0][1:]] == sorted(PypeItDataPaths.defined_paths)


def test_detector_table(monkeypatch, tmp_path):
    script = load_script(monkeypatch, 'build_detector_table')
    output = tmp_path / 'inst_detector_table.rst'
    script.write_detector_table(output)
    rows = table_rows(output.read_text())
    assert len(rows) == 1
    assert rows[0][0][:3] == ['Instrument', 'Det', 'specaxis']
    assert any(row[0] == 'keck_lris_blue' for row in rows[0][1:])


def test_image_processing_defaults_table(monkeypatch, tmp_path):
    script = load_script(monkeypatch, 'build_imgproc_defaults')
    output = tmp_path / 'imgproc_defaults_table.rst'
    script.write_imgproc_def_table(output)
    rows = table_rows(output.read_text())
    assert len(rows) == 1
    assert rows[0][0][:2] == ['Parameter', 'Default']
    assert 'apply_gain' in {row[0] for row in rows[0][1:]}


def test_spectrograph_table(monkeypatch, tmp_path):
    script = load_script(monkeypatch, 'build_spectbl_rst')
    script.write_spec_table(tmp_path)
    rows = table_rows((tmp_path / 'spectrographs_table.rst').read_text())
    assert len(rows) == 1
    assert rows[0][0][:3] == ['PypeIt Name', 'PypeIt Class', 'Telescope']
    assert 'keck_lris_blue' in {row[0] for row in rows[0][1:]}


def test_standard_star_tables(monkeypatch, tmp_path):
    script = load_script(monkeypatch, 'build_standardstbl_rst')
    from pypeit.core import standard

    script.write_tables(tmp_path)
    for archive in [*standard.get_archive_sets(), 'blackbody', 'kurucz93']:
        rows = table_rows((tmp_path / f'{archive}_table.rst').read_text())
        assert len(rows) == 1
        assert len(rows[0]) > 1


def test_telluric_grid_table(monkeypatch, tmp_path):
    script = load_script(monkeypatch, 'build_telluric_table')
    monkeypatch.setattr(script.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(
        returncode=0, stdout='2026-01-01 00:00:00 2.0 GiB example.fits\n', stderr=''))
    monkeypatch.setitem(script.os.environ, 'TZ', script.os.environ.get('TZ', 'UTC'))

    output = tmp_path / 'TelFit_files.rst'
    script.tellgrid_table('TelFit', output)
    rows = table_rows(output.read_text())
    assert rows == [[['File', 'Size', 'Last Modified (UTC)'],
                     ['example.fits', '2.0 GiB', '2026-01-01 00:00:00']]]
