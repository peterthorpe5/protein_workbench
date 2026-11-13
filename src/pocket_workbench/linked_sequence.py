"""Validate coordinates for the linked Nightingale and molecular viewer."""

from __future__ import annotations

import hashlib
import re
from functools import lru_cache
from importlib.resources import files

MAX_AXIS_LENGTH = 100_000
FASTA_MAPPING = {"exact_explicit", "exact_sequence", "unique_exact_fragment_explicit"}


def mapped_positions(
    *, residues: list[dict], sequence: str, field: str
) -> list[int] | None:
    """Check a contiguous exact sequence mapping to a given coordinate field.

    Args:
        residues: Prepared model's ordered residue records.
        sequence: One-based axis sequence, with position one at its first letter.
        field: Coordinate key used for this axis.

    Returns:
        Verified consecutive axis positions, or None if any row disagrees.
    """
    if not residues or not sequence or len(sequence) > MAX_AXIS_LENGTH:
        return None
    start = residues[0].get(field)
    if type(start) is not int or not 1 <= start <= len(sequence) - len(residues) + 1:
        return None
    positions = list(range(start, start + len(residues)))
    for position, row in zip(positions, residues, strict=True):
        if row.get(field) != position or type(row.get(field)) is not int:
            return None
        if row.get("amino_acid") != sequence[position - 1]:
            return None
    return positions


@lru_cache(maxsize=1)
def recorded_landmarks() -> dict[str, tuple[str, list[dict]]]:
    """Load the BRCA1 project review windows keyed by exact sequence hash.

    Returns:
        A sequence hash to (accession, labelled feature ranges) lookup.
    """
    data = (
        files("pocket_workbench")
        .joinpath("content", "BRCA1_LINKED_LANDMARKS.tsv")
        .read_text(encoding="utf-8")
    )
    landmarks: dict[str, tuple[str, list[dict]]] = {}
    for line in data.splitlines()[1:]:
        digest, accession, start, end, label, kind = line.split("\t")
        _, features = landmarks.setdefault(digest, (accession, []))
        features.append(
            {"start": int(start), "end": int(end), "label": label, "kind": kind}
        )
    return landmarks


def linked_sequence_axis(
    *, residues: list[dict], model: dict, representative_sequence: str | None = None
) -> dict | None:
    """Use verified FASTA positions, otherwise fall back to prepared-model positions.

    Args:
        residues: Existing portable residue map, in prepared-model order.
        model: Model metadata, including its recorded FASTA mapping status.
        representative_sequence: Exact imported FASTA entry, when available.

    Returns:
        Axis sequence, coordinate field and sequence-verified landmarks, or None
        if the resource lacks a usable residue map.
    """
    model_sequence = "".join(str(row.get("amino_acid", "")) for row in residues)
    if not re.fullmatch(r"[A-Z]+", model_sequence) or len(model_sequence) != len(
        residues
    ):
        return None
    if not mapped_positions(
        residues=residues, sequence=model_sequence, field="model_position"
    ):
        return None
    keys = [row.get("residue_key") for row in residues]
    tool_positions = [row.get("tool_residue_number") for row in residues]
    if (
        any(not isinstance(key, str) or not key for key in keys)
        or len(set(keys)) != len(keys)
        or any(type(position) is not int or position < 1 for position in tool_positions)
        or len(set(tool_positions)) != len(tool_positions)
    ):
        return None
    sequence = model_sequence
    field = "model_position"
    label = "Prepared model positions"
    features: list[dict] = []
    if (
        representative_sequence
        and len(representative_sequence) <= MAX_AXIS_LENGTH
        and re.fullmatch(r"[A-Z]+", representative_sequence)
        and model.get("mapping_status") in FASTA_MAPPING
    ):
        positions = mapped_positions(
            residues=residues,
            sequence=representative_sequence,
            field="sequence_position",
        )
        if positions and (
            model.get("mapping_status") == "unique_exact_fragment_explicit"
            or model_sequence == representative_sequence
        ):
            sequence = representative_sequence
            field = "sequence_position"
            label = "Representative FASTA positions"
    digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    if digest in recorded_landmarks():
        accession, recorded = recorded_landmarks()[digest]
        features = [dict(feature) for feature in recorded]
        if field == "sequence_position":
            label += f" · {accession}"
    return {"sequence": sequence, "field": field, "label": label, "features": features}
