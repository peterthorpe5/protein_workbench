"""Check public reference and species UI behaviour on a native resource."""

from pathlib import Path
from unittest.mock import Mock

from streamlit.testing.v1 import AppTest

from pocket_workbench import discovery_ui, ui
from pocket_workbench.core import WorkbenchError, load_project, write_json


def test_public_reference_lookup_is_explicit_and_clears_stale_results(
    project, monkeypatch, tmp_path
):
    """A coordinate-free model offers lookup without changing its saved evidence."""
    data = load_project(root=project)
    model = data["models"][0]
    model.update(prepared_model=None, coordinates=None, source_accession="NP_001013434")
    path = project / "project.json"
    write_json(path=path, value=data)
    before = path.read_bytes()
    seq = data["sequences"][model["sequence_id"]]
    fasta = f">NP_001013434.1 [Canis lupus familiaris]\n{seq}\n".encode()
    fetch = Mock(return_value=({"accession": "NP_001013434.1", "sequence": seq}, fasta))
    afdb = Mock()
    monkeypatch.setattr(discovery_ui, "fetch_ncbi_fasta", fetch)
    monkeypatch.setattr(discovery_ui, "fetch_afdb_model", afdb)
    monkeypatch.setenv("POCKET_WORKBENCH_RESOURCE", str(project))
    monkeypatch.setenv("POCKET_WORKBENCH_WORKSPACE", str(tmp_path / "workspace"))
    app = AppTest.from_file(str(Path(ui.__file__)), default_timeout=20).run()
    assert not app.exception and not app.error
    fetch.assert_not_called()
    afdb.assert_not_called()
    app.button(key="model_reference_ncbi_fetch").click().run()
    fetch.assert_called_once_with(accession="NP_001013434")
    assert app.session_state["model_reference_ncbi_download"]["payload"] == fasta
    input_box = next(x for x in app.text_input if x.label == "NCBI protein accession")
    input_box.set_value("NP_001394510.1").run()
    assert "model_reference_ncbi_download" not in app.session_state
    fetch.side_effect = WorkbenchError("No public record was found")
    app.button(key="model_reference_ncbi_fetch").click().run()
    assert any("No public record" in row.value for row in app.warning)
    assert "model_reference_ncbi_download" not in app.session_state
    assert path.read_bytes() == before


def test_alignment_composition_and_clade_chart_use_the_saved_members(
    project, monkeypatch, tmp_path
):
    """The summary and selected clade use representative counts, with downloads."""
    data = load_project(root=project)
    for model in data["models"]:
        model["species"] = (
            "Fixture one"
            if model["sample_id"] == "synthetic_reference"
            else "Fixture two"
        )
    write_json(path=project / "project.json", value=data)
    monkeypatch.setenv("POCKET_WORKBENCH_RESOURCE", str(project))
    monkeypatch.setenv("POCKET_WORKBENCH_WORKSPACE", str(tmp_path / "workspace"))
    app = AppTest.from_file(str(Path(ui.__file__)), default_timeout=20).run()
    app.sidebar.radio[0].set_value("Alignment & trees").run()
    assert not app.exception and not app.error
    assert any(
        "Download species counts" in row.proto.label
        for row in app.get("download_button")
    )
    app.checkbox(key="saved_alignment_tree_preview").check().run()
    assert not app.exception and not app.error
    assert app.selectbox(key="saved_alignment_subcluster_compare").options
    assert len(app.get("plotly_chart")) >= 2
    assert any(
        "Download subcluster species counts" in row.proto.label
        for row in app.get("download_button")
    )
