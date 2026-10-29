"""Check the AstroPy-formatted tables produced by the documentation scripts."""

import importlib
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

from docutils import nodes
from docutils.core import publish_doctree
import numpy as np
import pytest


DOC_SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
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


def test_astropy_rst_writer(monkeypatch):
    writer = load_script(monkeypatch, 'rst_table')
    data = np.array([
        ['Name', 'Description', 'Optional'],
        ['``first``', 'A :math:`x` value', ''],
        ['``second``', 'Another value', 'present'],
    ], dtype=object)

    source = writer.format_rst_table(data)
    rows = table_rows(source)
    assert len(rows) == 1
    assert rows[0][0] == ['Name', 'Description', 'Optional']
    assert [row[0] for row in rows[0][1:]] == ['first', 'second']
    assert rows[0][1][2] == ''
    assert ':math:`x`' in source
    assert all(line == line.rstrip() for line in source.splitlines())


def test_headerless_rst_writer(monkeypatch):
    writer = load_script(monkeypatch, 'rst_table')
    data = np.array([['First', '``one``'], ['Second', '``two``']], dtype=object)

    source = writer.format_rst_table(data, has_header=False)
    assert table_rows(source) == [[['First', 'one'], ['Second', 'two']]]
    assert not list(publish_doctree(source).findall(nodes.thead))
    assert all(line == line.rstrip() for line in source.splitlines())


def test_afterburn_datamodel_tables(monkeypatch, tmp_path):
    script = load_script(monkeypatch, 'build_afterburn_datamodels_rst')
    script.sens_datamodel(tmp_path)
    script.telluric_datamodel(tmp_path)

    sens = table_rows((tmp_path / 'datamodel_sensfunc.rst').read_text())
    telluric = table_rows((tmp_path / 'datamodel_telluric.rst').read_text())
    assert len(sens) == 3
    assert len(telluric) == 2
    assert sens[0][0] == ['HDU Name', 'HDU Type', 'Data Type', 'Description']
    assert {'PRIMARY', 'TELLURIC', 'SENS', 'WAVE', 'ZEROPOINT', 'THROUGHPUT'} \
        <= {row[0] for row in sens[0][1:]}
    assert telluric[1][0] == ['Column', 'Data Type', 'Description']


def test_calibration_image_tables(monkeypatch, tmp_path):
    script = load_script(monkeypatch, 'build_calib_images_rst')
    from pypeit.images.buildimage import ArcImage

    arc = table_rows('\n'.join(script.basic_pypeitimage_datamodel(ArcImage)))
    assert len(arc) == 1
    assert arc[0][0] == ['HDU Name', 'HDU Type', 'Data Type', 'Description']
    assert {row[0] for row in arc[0][1:]} == {'PRIMARY', *ArcImage.output_to_disk}

    script.slits_datamodel(tmp_path)
    slits = table_rows((tmp_path / 'datamodel_slittraceset.rst').read_text())
    assert len(slits) == 2
    assert slits[1][0] == ['Column', 'Data Type', 'Description']
    assert 'TRACEID' in {row[0] for row in slits[1][1:]}


@pytest.mark.parametrize('class_name', ['SpecObj', 'Spec2DObj', 'OneSpec'])
def test_specobj_datamodel_tables(monkeypatch, tmp_path, class_name):
    script = load_script(monkeypatch, 'build_specobj_rst')
    modules = {'SpecObj': script.specobj, 'Spec2DObj': script.spec2dobj,
               'OneSpec': script.coadd1d}
    obj = getattr(modules[class_name], class_name)
    script.write_specobj_table(obj, tmp_path)

    rows = table_rows((tmp_path / f'datamodel_{class_name.lower()}.rst').read_text())
    assert len(rows) == 1
    assert rows[0][0] == ['Obj Key', 'Obj Type', 'Array Type', 'Description']
    assert {row[0] for row in rows[0][1:]} == set(obj.datamodel)


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
    from pypeit.pkg.pypeitdata import PypeItDataPaths

    output = tmp_path / 'data_dir.rst'
    script.write_table(output)
    rows = table_rows(output.read_text())
    assert len(rows) == 1
    assert rows[0][0] == ['Reference', 'Subdirectory', 'Host']
    assert [row[0] for row in rows[0][1:]] == sorted(PypeItDataPaths.defined_paths)


def test_datacontainer_datamodel_table(monkeypatch):
    script = load_script(monkeypatch, 'build_datacontainer_datamodels')
    from pypeit.images.detector_container import DetectorContainer

    source = script.build_datamodel_tbl(DetectorContainer)
    rows = table_rows(source)
    assert len(rows) == 1
    assert rows[0][0] == ['Attribute', 'Type', 'Array Type', 'Description']
    assert {row[0] for row in rows[0][1:]} == set(DetectorContainer.datamodel)
    assert all(line == line.rstrip() for line in source.splitlines())


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


def test_contributors_table(monkeypatch, tmp_path):
    script = load_script(monkeypatch, 'write_contributors_table')
    doc_root = tmp_path / 'doc'
    (doc_root / 'scripts').mkdir(parents=True)
    (doc_root / 'include').mkdir()
    (doc_root / 'scripts' / 'contributors.csv').write_text(
        'first,First Contributor,True\nsecond,Second Contributor,False\n')
    monkeypatch.setattr(script.resources, 'files', lambda _: tmp_path / 'pypeit')

    script.main()

    rows = table_rows((doc_root / 'include' / 'contributors.rst').read_text())
    assert rows == [[['Name', 'GitHub handle'], ['First Contributor', '@first']]]


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


def test_dependency_table(monkeypatch, tmp_path):
    script = load_script(monkeypatch, 'build_dependency_rst')
    script.write_dependency_table(PROJECT_FILE, tmp_path)
    source = (tmp_path / 'dependencies_table.rst').read_text()
    rows = table_rows(source)

    assert len(rows) == 1
    assert [row[0] for row in rows[0]] == [
        'Python Version', 'Required for users', 'Required for developers']
    assert 'astropy' in rows[0][1][1]
    assert 'pytest' in rows[0][2][1]
    assert not list(publish_doctree(source).findall(nodes.thead))
    assert all(line == line.rstrip() for line in source.splitlines())


def test_parset_rst_tables():
    from pypeit.par.parset import ParSet

    class ChildPar(ParSet):
        pass

    child = ChildPar(['beta'], values=[2], dtypes=[int])
    par = ParSet(['process', 'alpha'], values=[child, 1],
                 dtypes=[ChildPar, int], descr=['Image processing', 'A value'])
    source = '\n'.join(par.to_rst_table())
    rows = table_rows(source)

    assert len(rows) == 2
    assert rows[0][0] == ['Key', 'Type', 'Options', 'Default', 'Description']
    assert {row[0] for row in rows[0][1:]} == {'alpha', 'process'}
    assert rows[1][0] == ['Key', 'Type', 'Options', 'Default', 'Description']
    assert {row[0] for row in rows[1][1:]} == set(child.keys())
