"""Test exact public IDs, allowlisted downloads and saved sequence checks offline."""

from __future__ import annotations

import hashlib
import json
from urllib.error import HTTPError

import gemmi
import pytest

from pocket_workbench import public_records as records
from pocket_workbench.core import WorkbenchError


def test_ncbi_accepts_full_protein_accessions_and_checks_the_returned_version(
    monkeypatch,
):
    """Resolve unversioned lookups while keeping versioned lookups exact."""
    assert records.ncbi_protein_accession(value=" np_001013434.1 ") == "NP_001013434.1"
    assert records.ncbi_protein_accession(value="KAG3268860.1") == "KAG3268860.1"
    for value in (
        None,
        "XM_123456",
        "fold_np_001013434_model_0",
        "P12345",
        "1&db=nuccore",
    ):
        assert records.ncbi_protein_accession(value=value) is None
    payload = b">NP_001013434.1 BRCA1 [Canis lupus familiaris]\nACD\nEF\n"
    calls = []

    def supplied(*, url, hosts, maximum):
        calls.append((url, hosts, maximum))
        return payload

    monkeypatch.setattr(records, "download", supplied)
    result, actual = records.fetch_ncbi_fasta(accession="NP_001013434")
    assert (
        result["species"] == "Canis lupus familiaris" and result["sequence"] == "ACDEF"
    )
    assert result["accession"] == "NP_001013434.1" and actual == payload
    assert result["sha256"] == hashlib.sha256(payload).hexdigest()
    assert "db=protein&id=NP_001013434&rettype=fasta&retmode=text" in calls[0][0]
    for requested in ("NP_001013434.2", "NP_001394510"):
        with pytest.raises(WorkbenchError, match="different accession"):
            records.parse_ncbi_fasta(payload=payload, requested=requested)
    with pytest.raises(WorkbenchError, match="exactly one"):
        records.parse_ncbi_fasta(payload=payload * 2, requested="NP_001013434")


@pytest.mark.parametrize(
    "url",
    [
        "http://example.org/a",
        "https://example.org.evil.test/a",
        "https://a@example.org/a",
        "https://example.org:8080/a",
        "https://example.org/a#b",
        "https://example.org/a b",
        "https://example.org\\a",
    ],
)
def test_external_url_and_redirect_reject_unapproved_destinations(url):
    """Metadata responses and redirects cannot select a different service."""
    hosts = frozenset({"example.org"})
    with pytest.raises(WorkbenchError, match="unsupported"):
        records.validated_url(url=url, hosts=hosts)
    with pytest.raises(WorkbenchError, match="unsupported"):
        records.RestrictedRedirect(hosts=hosts).redirect_request(
            None, None, 302, "", {}, url
        )


def test_download_bounds_and_http_errors(monkeypatch):
    """Transport failures are explicit and oversized responses are not retained."""
    url = "https://example.org/protein"
    opts = {"url": url, "hosts": frozenset({"example.org"}), "maximum": 3}

    class Reply:
        headers = {}

        def __enter__(self):
            return self

        def __exit__(self, *_):
            return None

        def read(self, maximum):
            return b"abcd"[:maximum]

    class Opener:
        def open(self, request, timeout):
            return Reply()

    monkeypatch.setattr(records, "build_opener", lambda *_: Opener())
    with pytest.raises(WorkbenchError, match="exceeded"):
        records.download(**opts)

    def unavailable(self, request, timeout):
        raise HTTPError(url, 404, "Missing", {}, None)

    monkeypatch.setattr(Opener, "open", unavailable)
    with pytest.raises(WorkbenchError, match="No public record"):
        records.download(**opts)


def test_afdb_binds_model_metadata_and_atomic_sequence(project, monkeypatch):
    """Return a current-schema complete model; reject ambiguity or mismatched atoms."""
    pdb = project / "models" / "synthetic_reference" / "prepared.pdb"
    structure = gemmi.read_structure(str(pdb))
    seq = "".join(
        gemmi.find_tabulated_residue(r.name).one_letter_code for r in structure[0][0]
    )
    payload = structure.make_mmcif_document().as_string().encode()
    metadata = {
        "modelEntityId": "AF-P12345-F1",
        "uniprotAccession": "P12345",
        "sequenceStart": 1,
        "sequenceEnd": len(seq),
        "sequence": seq,
        "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v6.cif",
    }
    calls = []

    def make(*rows):
        """Encode the mocked AFDB list response."""
        return json.dumps(rows).encode()

    def supplied(*, url, hosts, maximum):
        calls.append(url)
        return make(metadata) if "/api/" in url else payload

    monkeypatch.setattr(records, "download", supplied)
    record, body = records.fetch_afdb_model(accession="P12345", expected=seq)
    assert record["model_id"] == "AF-P12345-F1" and record["sequence"] == seq
    assert body == payload and len(calls) == 2
    viewer = records.public_model_document(record=record, payload=body)
    assert "Public AlphaFold DB sequence positions" in viewer
    assert '"pockets": []' in viewer and "synthetic_reference" not in viewer
    with pytest.raises(WorkbenchError, match="differ from their retrieval record"):
        records.public_model_document(record=record, payload=body + b"\n")
    with pytest.raises(WorkbenchError, match="No unique"):
        records.fetch_afdb_model(accession="P12345", expected="AAAA")
    with pytest.raises(WorkbenchError, match="No unique"):
        records.afdb_choice(
            payload=make(metadata, metadata), accession="P12345", expected=seq
        )
    with pytest.raises(WorkbenchError, match="does not match"):
        records.validate_afdb_coordinates(payload=payload, sequence="AAAA")
