"""Coordinate validation and offline Nightingale integration tests."""

from __future__ import annotations

import hashlib
import os

import pytest

from pocket_workbench.core import load_project
from pocket_workbench.linked_sequence import (
    linked_sequence_axis,
    mapped_positions,
    recorded_landmarks,
)
from pocket_workbench.viewers import (
    model_viewer_records,
    structure_document,
    viewer_asset,
)


def residue_map(*, sequence: str, offset: int = 0) -> list[dict]:
    """Provide a small prepared numbering map with distinct parent positions."""
    return [
        {
            "residue_key": f"A_{position}",
            "tool_residue_number": position,
            "model_position": position,
            "sequence_position": position + offset,
            "amino_acid": aa,
        }
        for position, aa in enumerate(sequence, 1)
    ]


def test_exact_fragment_maps_to_parent_and_invalid_rows_fall_back():
    """Verify parent coordinates by sequence and require an exact mapping status."""
    rows = residue_map(sequence="CDE", offset=2)
    model = {"mapping_status": "unique_exact_fragment_explicit"}
    sequence = "AACDEF"
    assert mapped_positions(
        residues=rows, sequence=sequence, field="sequence_position"
    ) == [3, 4, 5]
    axis = linked_sequence_axis(
        residues=rows, model=model, representative_sequence=sequence
    )
    assert axis["field"] == "sequence_position" and axis["sequence"] == sequence
    ambiguous = linked_sequence_axis(
        residues=rows,
        model={"mapping_status": "no_exact_fasta_match"},
        representative_sequence=sequence,
    )
    assert ambiguous["field"] == "model_position" and ambiguous["sequence"] == "CDE"
    rows[1]["sequence_position"] = 5
    invalid = linked_sequence_axis(
        residues=rows, model=model, representative_sequence=sequence
    )
    assert invalid["field"] == "model_position"
    rows[1]["residue_key"] = rows[0]["residue_key"]
    assert linked_sequence_axis(residues=rows, model=model) is None


def test_recorded_landmarks_require_an_exact_sequence_hash(monkeypatch):
    """Ensure a named protein or same-length substitution cannot add landmarks."""
    known = recorded_landmarks()
    dog = next(
        features
        for accession, features in known.values()
        if accession == "NP_001013434"
    )
    assert [(x["start"], x["end"]) for x in dog] == [
        (1653, 1727),
        (1757, 1761),
        (1762, 1766),
        (1767, 1850),
    ]
    seq = "ACDEF"
    digest = hashlib.sha256(seq.encode()).hexdigest()
    landmark = {"start": 2, "end": 3, "kind": "domain", "label": "Fixture domain"}
    monkeypatch.setattr(
        "pocket_workbench.linked_sequence.recorded_landmarks",
        lambda: {digest: ("test", [landmark])},
    )
    exact = linked_sequence_axis(residues=residue_map(sequence=seq), model={})
    other = linked_sequence_axis(residues=residue_map(sequence="ACGEF"), model={})
    assert exact["features"] == [landmark]
    assert other["features"] == []


def test_linked_assets_are_bundled_and_not_added_to_pair_views(project):
    """Single-model downloads include Nightingale without changing pair overlays."""
    data = load_project(root=project)
    model = data["models"][0]
    metadata = model_viewer_records(
        root=project, model=model, pockets=data["pockets"], sequences=data["sequences"]
    )
    assert metadata["linked_axis"]["field"] == "sequence_position"
    pdb = (project / model["prepared_model"]).read_text()
    document = structure_document(
        reference_pdb=pdb, reference_label="test", reference_metadata=metadata
    )
    assert viewer_asset(name="nightingale-components.js").strip()[:80] in document
    assert 'src="http' not in document
    paired = structure_document(
        reference_pdb=pdb,
        reference_label="test",
        reference_metadata=metadata,
        mobile_pdb=pdb,
        mobile_label="test mobile",
        mobile_metadata=metadata,
    )
    assert viewer_asset(name="nightingale-components.js").strip()[:80] not in paired


@pytest.mark.skipif(
    os.getenv("POCKET_WORKBENCH_BROWSER_TESTS") != "1",
    reason="Requires Playwright Chromium",
)
def test_offline_nightingale_selection_reaches_the_molecular_viewer(project, tmp_path):
    """Exercise real pocket/sequence clicks, missing positions and display colours."""
    from playwright.sync_api import sync_playwright

    data = load_project(root=project)
    model = next(
        row for row in data["models"] if row["sample_id"] == "synthetic_reference"
    )
    metadata = model_viewer_records(
        root=project, model=model, pockets=data["pockets"], sequences=data["sequences"]
    )
    metadata["linked_axis"]["sequence"] = "E" * 6 + metadata["linked_axis"]["sequence"]
    metadata["linked_axis"]["features"] = [
        {"start": 7, "end": 11, "label": "Fixture domain", "kind": "domain"}
    ]
    for row in metadata["residues"]:
        row["sequence_position"] += 6
    path = tmp_path / "linked.html"
    path.write_text(
        structure_document(
            reference_pdb=(project / model["prepared_model"]).read_text(),
            reference_label=model["sample_id"],
            reference_metadata=metadata,
        )
    )
    with sync_playwright() as api:
        browser = api.chromium.launch(
            args=[
                "--no-sandbox",
                "--enable-unsafe-swiftshader",
                "--use-angle=swiftshader",
            ]
        )
        page = browser.new_page(viewport={"width": 1400, "height": 1100})
        errors, network = [], []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.on(
            "request",
            lambda request: (
                network.append(request.url) if request.url.startswith("http") else None
            ),
        )
        page.route("http**://**", lambda route: route.abort())
        page.goto(path.as_uri())
        page.wait_for_function("Boolean(window.pocketDisplayState)")
        assert page.locator("#linked-sequence").is_visible()
        assert (
            "Representative FASTA positions"
            in page.locator("#linked-sequence").inner_text()
        )
        search = page.get_by_role(
            "searchbox", name="Filter pocket tracks by native ID or rank"
        )
        search.fill("2")
        row = page.locator(".sequence-pocket-rows .sequence-track-row").first
        row.locator("nightingale-track rect.outer-rectangle").first.click()
        assert page.evaluate("window.pocketDisplayState[0].pockets") == [2]
        page.get_by_role("button", name="Fixture domain").click()
        assert (
            page.evaluate("window.pocketDisplayState[0].landmark") == "Fixture domain"
        )
        position = page.get_by_label("Residue position", exact=True)
        position.fill("12")
        page.get_by_role("button", name="Find residue").click()
        assert page.evaluate("window.pocketDisplayState[0].picked") == 6
        assert "Original A:6" in page.locator("#linked-selection").inner_text()
        assert "FASTA position: 12" in page.locator("#picked").inner_text()
        assert (
            page.locator("nightingale-sequence").get_attribute("highlight") == "12:12"
        )
        seq = page.locator("nightingale-sequence text.base").filter(has_text="H").first
        point = seq.bounding_box()
        assert point is not None
        page.mouse.click(
            point["x"] + point["width"] / 2, point["y"] + point["height"] / 2
        )
        assert page.evaluate("window.pocketDisplayState[0].picked") == 7
        page.evaluate("""() => {
          const model = window.pocketViewer.getModel(0);
          const [atom] = model.selectedAtoms({atom:'CA',resi:8});
          atom.callback(atom);
        }""")
        assert position.input_value() == "14"
        position.fill("2")
        page.get_by_role("button", name="Find residue").click()
        assert page.evaluate("window.pocketDisplayState[0].picked") is None
        assert "No coordinates" in page.locator("#linked-selection").inner_text()
        page.get_by_text("Pocket appearance", exact=True).click()
        page.get_by_label("Pocket colours identify").select_option("custom")
        page.get_by_label("Pocket colour", exact=True).fill("#5a3bc4")
        page.get_by_label("Pocket colour", exact=True).dispatch_event("change")
        assert page.evaluate("window.pocketDisplayState[0].colours") == ["#5a3bc4"]
        page.get_by_text("Domain and motif colours", exact=True).click()
        page.get_by_label("Domain colour", exact=True).fill("#385c17")
        page.get_by_label("Domain colour", exact=True).dispatch_event("change")
        assert (
            page.get_by_role("button", name="Fixture domain").evaluate(
                "e => e.style.getPropertyValue('--feature-colour')"
            )
            == "#385c17"
        )
        page.get_by_role("button", name="Clear highlights").click()
        state = page.evaluate("window.pocketDisplayState[0]")
        assert state["pockets"] == [2] and state["landmark"] is None
        assert not errors, errors
        assert not network, network
        browser.close()
