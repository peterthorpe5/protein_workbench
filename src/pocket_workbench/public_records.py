"""Fetch explicit public protein records without modifying saved evidence."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import re
import threading
import time
from datetime import UTC, datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

import gemmi

from .core import WorkbenchError
from .external_tools import canonical_uniprot

LOGGER = logging.getLogger(__name__)
AFDB = "alphafold.ebi.ac.uk"
NCBI = "eutils.ncbi.nlm.nih.gov"
MAX_FASTA = 250_000
MAX_METADATA = 2_000_000
MAX_CIF = 24_000_000
PROTEIN_ID = re.compile(
    r"(?:(?:AP|NP|XP|YP|WP|ZP)_[0-9]{5,12}|[A-Z]{3}[0-9]{5,12})(?:\.[0-9]{1,4})?"
)
FASTA_ALPHABET = re.compile(r"[ACDEFGHIKLMNPQRSTVWYBXZJUO]{1,100000}")
_NCBI_LOCK = threading.Lock()
_NCBI_NEXT = 0.0


def ncbi_protein_accession(*, value: object) -> str | None:
    """Accept only a whole supported RefSeq or GenBank protein accession."""
    if not isinstance(value, str):
        return None
    value = value.strip().upper()
    if PROTEIN_ID.fullmatch(value):
        return value
    return None


def validated_url(*, url: str, hosts: frozenset[str]) -> str:
    """Allow only HTTPS URLs on explicitly approved hosts and port 443."""
    try:
        if not isinstance(url, str):
            raise ValueError("invalid destination")
        parts = urlsplit(url)
        if (
            parts.scheme != "https"
            or parts.hostname not in hosts
            or parts.username is not None
            or parts.password is not None
            or parts.port not in {443, None}
            or parts.fragment
            or re.search(r"[\x00-\x20\\]", url)
        ):
            raise ValueError("invalid destination")
    except (ValueError, TypeError) as error:
        raise WorkbenchError(
            "The database returned an unsupported download URL"
        ) from error
    return url


class RestrictedRedirect(HTTPRedirectHandler):
    """Check every server redirect before urllib follows it."""

    def __init__(self, *, hosts: frozenset[str]) -> None:
        """Retain the allowlist for this request."""
        self.hosts = hosts

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        """Reject redirects that leave the approved service endpoints."""
        validated_url(url=newurl, hosts=self.hosts)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def download(*, url: str, hosts: frozenset[str], maximum: int) -> bytes:
    """Read a bounded public response using only HTTPS and approved redirects."""
    validated_url(url=url, hosts=hosts)
    opener = build_opener(RestrictedRedirect(hosts=hosts))
    request = Request(
        url,
        headers={"User-Agent": "ProteinPocketWorkbench/1.6 (public protein lookup)"},
    )
    try:
        with opener.open(request, timeout=12) as handle:
            declared = handle.headers.get("Content-Length")
            if declared and (not declared.isdecimal() or int(declared) > maximum):
                raise WorkbenchError(
                    "Database response exceeds the download size limit"
                )
            body = handle.read(maximum + 1)
    except HTTPError as error:
        if error.code == 404:
            raise WorkbenchError(
                "No public record was found for this exact accession"
            ) from error
        raise WorkbenchError(f"Database request returned HTTP {error.code}") from error
    except (URLError, TimeoutError, OSError) as error:
        LOGGER.warning("Public protein download failed: %s", type(error).__name__)
        raise WorkbenchError(
            "Could not reach the database. Try again or use its website."
        ) from error
    if not body or len(body) > maximum:
        raise WorkbenchError(
            "Database response was empty or exceeded the download size limit"
        )
    return body


def parse_ncbi_fasta(*, payload: bytes, requested: str) -> dict:
    """Validate exactly one NCBI FASTA record with the requested accession/version."""
    try:
        header, *body = payload.decode("utf-8").strip().splitlines()
    except (ValueError, UnicodeError) as error:
        raise WorkbenchError("NCBI did not return protein FASTA") from error
    returned = (header[1:].split() or [""])[0] if header.startswith(">") else ""
    seq = "".join(line.strip() for line in body).upper()
    if not ncbi_protein_accession(value=returned) or not FASTA_ALPHABET.fullmatch(seq):
        raise WorkbenchError("NCBI did not return exactly one valid protein FASTA")
    if requested != (returned if "." in requested else returned.split(".", 1)[0]):
        raise WorkbenchError("NCBI returned a different accession or sequence version")
    species = re.search(r"\[([^\[\]]{1,150})\]\s*$", header)
    return {
        "accession": returned,
        "sequence": seq,
        "species": species[1] if species else "",
    }


def fetch_ncbi_fasta(*, accession: str) -> tuple[dict, bytes]:
    """Retrieve a single NCBI protein FASTA and return its audited versioned record."""
    global _NCBI_NEXT
    accession = ncbi_protein_accession(value=accession)
    if not accession:
        raise WorkbenchError("Use a complete RefSeq or GenBank protein accession")
    params = {
        "db": "protein",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "tool": "ProteinPocketWorkbench",
    }
    contact = os.getenv("POCKET_WORKBENCH_NCBI_EMAIL", "").strip()
    if contact and re.fullmatch(r"[^\s@]{1,100}@[^\s@]{1,100}", contact):
        params["email"] = contact
    query = urlencode(params)
    url = f"https://{NCBI}/entrez/eutils/efetch.fcgi?{query}"
    with _NCBI_LOCK:
        time.sleep(max(0.0, _NCBI_NEXT - time.monotonic()))
        _NCBI_NEXT = time.monotonic() + 0.4
    payload = download(url=url, hosts=frozenset({NCBI}), maximum=MAX_FASTA)
    record = parse_ncbi_fasta(payload=payload, requested=accession)
    record.update(
        source_url=f"https://www.ncbi.nlm.nih.gov/protein/{record['accession']}",
        sha256=hashlib.sha256(payload).hexdigest(),
        retrieved_utc=datetime.now(UTC).isoformat(),
    )
    return record, payload


def afdb_choice(*, payload: bytes, accession: str, expected: str | None) -> dict:
    """Choose one exact full-length AlphaFold DB prediction."""
    try:
        rows = json.loads(payload)
    except (ValueError, UnicodeError) as error:
        raise WorkbenchError("AlphaFold DB did not return valid metadata") from error
    if not isinstance(rows, list):
        raise WorkbenchError("AlphaFold DB did not return a prediction list")
    choices = []
    for row in rows:
        if not isinstance(row, dict) or row.get("uniprotAccession") != accession:
            continue
        sequence = row.get("sequence", row.get("uniprotSequence"))
        start = row.get("sequenceStart", row.get("uniprotStart"))
        end = row.get("sequenceEnd", row.get("uniprotEnd"))
        model_id = row.get("modelEntityId", row.get("entryId"))
        if (
            not isinstance(sequence, str)
            or not FASTA_ALPHABET.fullmatch(sequence)
            or start != 1
            or end != len(sequence)
            or (expected and expected != sequence)
            or not isinstance(model_id, str)
            or not re.fullmatch(r"[A-Za-z0-9_.-]{1,100}", model_id)
            or not isinstance(row.get("cifUrl"), str)
        ):
            continue
        choices.append(
            dict(
                accession=accession,
                sequence=sequence,
                model_id=model_id,
                source_url=validated_url(
                    url=row["cifUrl"], hosts=frozenset({AFDB, "alphafold.com"})
                ),
            )
        )
    if len(choices) != 1:
        raise WorkbenchError(
            "No unique full-length AlphaFold DB model matches the selected sequence"
        )
    return choices[0]


def validate_afdb_coordinates(*, payload: bytes, sequence: str) -> None:
    """Verify a single protein chain agrees with the declared sequence."""
    try:
        doc = gemmi.cif.read_string(payload.decode("utf-8"))
        structure = gemmi.make_structure_from_block(doc.sole_block())
        if len(structure) != 1 or len(structure[0]) != 1:
            raise ValueError("Expected one model and one chain")
        residues = [
            r
            for r in structure[0][0]
            if gemmi.find_tabulated_residue(r.name).is_amino_acid()
        ]
        observed = "".join(
            gemmi.find_tabulated_residue(r.name).one_letter_code for r in residues
        ).upper()
        atoms = [residue.find_atom("CA", "*") for residue in residues]
        if sequence != observed or any(
            not atom
            or not all(math.isfinite(value) for value in atom.pos)
            or not 0 <= atom.b_iso <= 100
            for atom in atoms
        ):
            raise ValueError("Coordinates or sequence differ")
    except (RuntimeError, ValueError, UnicodeError) as error:
        raise WorkbenchError(
            "The downloaded model does not match its declared sequence"
        ) from error


def fetch_afdb_model(
    *, accession: str, expected: str | None = None
) -> tuple[dict, bytes]:
    """Retrieve and validate an exact public model as a separate mmCIF download."""
    accession = canonical_uniprot(value=accession)
    if not accession:
        raise WorkbenchError("AlphaFold DB lookup needs a canonical UniProt accession")
    allowed = frozenset({AFDB, "alphafold.com"})
    raw = download(
        url=f"https://{AFDB}/api/prediction/{accession}",
        hosts=allowed,
        maximum=MAX_METADATA,
    )
    record = afdb_choice(payload=raw, accession=accession, expected=expected)
    payload = download(url=record["source_url"], hosts=allowed, maximum=MAX_CIF)
    validate_afdb_coordinates(payload=payload, sequence=record["sequence"])
    record.update(
        sha256=hashlib.sha256(payload).hexdigest(),
        retrieved_utc=datetime.now(UTC).isoformat(),
    )
    return record, payload


def public_model_document(*, record: dict, payload: bytes) -> str:
    """Render a separate public model with its own verified axis and no pocket rows.

    Args:
        record: Validated AlphaFold DB record, including exact sequence and SHA-256.
        payload: Unchanged mmCIF bytes from that record.

    Returns:
        Standalone offline HTML for the public model, ready for motif search.
    """
    from .importing import normalise_model
    from .linked_sequence import linked_sequence_axis
    from .viewers import structure_document

    if hashlib.sha256(payload).hexdigest() != record["sha256"]:
        raise WorkbenchError("Public model bytes differ from their retrieval record")
    validate_afdb_coordinates(payload=payload, sequence=record["sequence"])
    with TemporaryDirectory(prefix="workbench-reference-") as directory:
        source = Path(directory) / "reference.cif"
        destination = source.with_suffix(".pdb")
        source.write_bytes(payload)
        rows, _ = normalise_model(
            source=source,
            destination=destination,
            sample_id=record["model_id"],
            chain="",
            confidence="plddt",
        )
        pdb = destination.read_text(encoding="utf-8")
    for index, row in enumerate(rows, 1):
        row.update(model_position=index, sequence_position=index)
    if "".join(row["amino_acid"] for row in rows) != record["sequence"]:
        raise WorkbenchError(
            "Prepared public model does not preserve its full sequence"
        )
    axis = linked_sequence_axis(
        residues=rows,
        model={"mapping_status": "exact_explicit"},
        representative_sequence=record["sequence"],
    )
    if axis is None or axis["field"] != "sequence_position":
        raise WorkbenchError("A verified public model sequence axis could not be built")
    axis["label"] = f"Public AlphaFold DB sequence positions · {record['accession']}"
    return structure_document(
        reference_pdb=pdb,
        reference_label=f"Public AlphaFold DB · {record['model_id']}",
        reference_metadata={
            "confidence": "plddt",
            "residues": rows,
            "pockets": [],
            "linked_axis": axis,
            "initial_display": "none",
        },
        summary="Public database model. This view has no pocket predictions.",
    )
