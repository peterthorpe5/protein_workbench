"""Check motif matching, structural overlays and the downloaded scene PDF."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from pocket_workbench.core import load_project
from pocket_workbench.viewers import model_viewer_records, structure_document

pytestmark = pytest.mark.skipif(
    os.getenv("POCKET_WORKBENCH_BROWSER_TESTS") != "1",
    reason="Requires Playwright Chromium",
)


@pytest.fixture
def offline_page():
    """Use real Chromium with WebGL and reject external requests or script errors."""
    from playwright.sync_api import sync_playwright

    with sync_playwright() as api:
        browser = api.chromium.launch(
            args=[
                "--no-sandbox",
                "--enable-unsafe-swiftshader",
                "--use-angle=swiftshader",
            ]
        )
        page = browser.new_page(
            viewport={"width": 1400, "height": 1050}, accept_downloads=True
        )
        errors, requests = [], []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.route(
            "http**://**",
            lambda route: (requests.append(route.request.url), route.abort()),
        )
        yield page
        browser.close()
        assert not errors, errors
        assert not requests, requests


def fixture_document(*, project: Path, pair: bool) -> str:
    """Build a synthetic overlay with repeated sequences and a partial model."""
    data = load_project(root=project)
    fit = json.loads((project / data["comparisons"][0]["result_path"]).read_text())
    by_id = {model["sample_id"]: model for model in data["models"]}
    meta = {}
    for name in (fit["reference"], fit["mobile"]):
        model = by_id[name]
        record = model_viewer_records(
            root=project,
            model=model,
            pockets=data["pockets"],
            sequences=data["sequences"],
        )
        axis = record["linked_axis"]
        axis["features"] = [
            {"start": 2, "end": 8, "kind": "domain", "label": "Fixture BRCT core"}
        ]
        if name == fit["reference"]:
            # A verified parent prefix and suffix can contain unmodelled search hits.
            axis["sequence"] = "KLMN" + axis["sequence"] + "KLMN"
            for row in record["residues"]:
                row["sequence_position"] += 4
        meta[name] = record
    reference = fit["reference"]
    mobile = fit["mobile"]
    return structure_document(
        reference_pdb=(project / by_id[reference]["prepared_model"]).read_text(),
        reference_label=reference,
        reference_metadata=meta[reference],
        mobile_pdb=(project / fit["mobile_aligned"]).read_text() if pair else "",
        mobile_label=mobile if pair else "",
        mobile_metadata=meta[mobile] if pair else None,
        match=fit["pocket_matches"][0] if pair else None,
    )


def test_motif_occurrences_coordinates_colours_and_pdf(project, tmp_path, offline_page):
    """Use a real viewer to select all and individual hits without clearing a pocket."""
    path = tmp_path / "motif.html"
    path.write_text(fixture_document(project=project, pair=False), encoding="utf-8")
    page = offline_page
    page.goto(path.as_uri())
    page.wait_for_function("Boolean(window.pocketDisplayState)")
    page.get_by_role("button", name="Fixture BRCT core").click()
    initial = page.evaluate("window.pocketDisplayState[0].pockets")
    query = page.get_by_label("Amino acid sequence")
    query.fill("kl mn")
    page.get_by_role("button", name="Find amino acids").click()
    assert "4 exact occurrences" in page.locator("#motif-result").inner_text()
    assert page.evaluate("window.pocketDisplayState[0].motif") == [
        9,
        10,
        11,
        12,
        29,
        30,
        31,
        32,
    ]
    occurrence = page.get_by_label("Amino acid occurrence")
    occurrence.select_option("0")
    assert "no coordinates" in page.locator("#motif-result").inner_text()
    assert page.evaluate("window.pocketDisplayState[0].motif") == []
    occurrence.select_option("1")
    picker = page.get_by_label("Motif search colour")
    picker.fill("#176829")
    picker.dispatch_event("change")
    assert page.evaluate("window.pocketDisplayState[0].motif") == [9, 10, 11, 12]
    assert page.evaluate("window.pocketDisplayState[0].pockets") == initial
    assert page.evaluate("window.pocketDisplayState[0].landmark") == "Fixture BRCT core"
    assert (
        page.locator("nightingale-track")
        .filter(has=page.locator("#g_motif-search"))
        .count()
        == 1
    )
    with page.expect_download() as requested:
        page.get_by_role("button", name="Download PDF", exact=True).click()
    output = tmp_path / "view.pdf"
    requested.value.save_as(output)
    content = output.read_bytes()
    assert content.startswith(b"%PDF-1.4") and content.endswith(b"%%EOF\n")
    assert b"/Subtype /Image" in content and b"KLMN" in content
    assert b"Fixture BRCT" in content and b"fpocket pocket" in content
    assert b"0.0902 0.4078 0.1608 rg" in content
    assert len(content) > 10_000
    page.get_by_role("button", name="Clear motif search").click()
    assert page.evaluate("window.pocketDisplayState[0].motif") == []
    assert page.evaluate("window.pocketDisplayState[0].pockets") == initial


def test_pair_motif_and_landmark_use_the_chosen_models(project, tmp_path, offline_page):
    """Search a saved overlay without transferring positions across proteins."""
    path = tmp_path / "pair.html"
    path.write_text(fixture_document(project=project, pair=True), encoding="utf-8")
    page = offline_page
    page.goto(path.as_uri())
    page.wait_for_function("Boolean(window.pocketDisplayState)")
    assert not page.locator("#linked-sequence").is_visible()
    assert page.evaluate(
        "window.pocketStructureFeatures.motifHits('AAAA', 'aaa').hits"
    ) == [{"start": 1, "end": 3}, {"start": 2, "end": 4}]
    page.get_by_label("synthetic_mobile recorded feature").select_option("0")
    page.get_by_label("Model for amino acid search").select_option("1")
    page.get_by_label("Amino acid sequence").fill("KLMN")
    page.get_by_role("button", name="Find amino acids").click()
    first, second = page.evaluate("window.pocketDisplayState")
    assert not first["motif"] and second["motif"] == [9, 10, 11, 12]
    assert second["landmark"] == "Fixture BRCT core"
    assert first["pockets"] and second["pockets"]
    page.get_by_label("Model for amino acid search").select_option("both")
    first, second = page.evaluate("window.pocketDisplayState")
    assert first["motif"] == [9, 10, 11, 12, 29, 30, 31, 32]
    assert second["motif"] == [9, 10, 11, 12]
    page.get_by_label("Amino acid sequence").fill("K.M")
    page.get_by_role("button", name="Find amino acids").click()
    assert "Search is exact" in page.locator("#motif-result").inner_text()
    assert not page.evaluate("window.pocketDisplayState[1].motif")
