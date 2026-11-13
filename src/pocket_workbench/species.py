"""Describe species representation without counting model replicates as sequences."""

from __future__ import annotations

import hashlib
import io
import re
from collections import Counter, defaultdict
from importlib.resources import files

from Bio import Phylo

from .core import WorkbenchError

UNKNOWN = "Unreported species"
CONFLICT = "Conflicting species labels"
REVIEWED_TAXA = {
    "NP_001013434": "Canis lupus familiaris (dog)",
    "NP_001394510": "Homo sapiens (human)",
}


def read_species_map(
    *, text: str, id_column: str = "sequence_id", label_column: str = "species"
) -> dict[str, str]:
    """Read a bounded two-column TSV with explicit sequence-level assignments.

    Args:
        text: TSV text with a header and two columns.
        id_column: Expected header for the first column.
        label_column: Expected header for the second column.

    Returns:
        Unique ID to species or subcluster label mapping.
    """
    if len(text) > 1_000_000:
        raise WorkbenchError("The assignment TSV exceeds the size limit")
    lines = [line for line in text.lstrip("\ufeff").splitlines() if line.strip()]
    if not lines or lines[0].split("\t") != [id_column, label_column]:
        raise WorkbenchError(f"Use a TSV headed {id_column} and {label_column}")
    mapping = {}
    for number, line in enumerate(lines[1:], 2):
        parts = [part.strip() for part in line.split("\t")]
        if (
            len(parts) != 2
            or not all(parts)
            or len(parts[0]) > 200
            or len(parts[1]) > 160
            or any(char.isspace() for char in parts[0])
            or any(ord(c) < 32 for c in parts[1])
            or parts[0] in mapping
            or parts[1] in {UNKNOWN, CONFLICT}
        ):
            raise WorkbenchError(f"Invalid or duplicate assignment on TSV row {number}")
        mapping[parts[0]] = parts[1]
    if not mapping:
        raise WorkbenchError("The assignment TSV has no data rows")
    return mapping


def species_assignments(
    *, data: dict, identifiers: list[str], supplied: dict[str, str]
) -> dict[str, str]:
    """Resolve species from explicit metadata or a reviewed accession and sequence.

    Args:
        data: Project models and representative sequences keyed by exact ID.
        identifiers: Representative IDs to describe.
        supplied: Optional explicit species labels keyed by representative ID.

    Returns:
        Species label, or an unknown or conflict label, for each requested ID.
    """
    wanted = set(identifiers)
    by_id: dict[str, set[str]] = defaultdict(set)
    accessions = {name: {re.sub(r"\.\d+$", "", name)} for name in wanted}
    for model in data["models"]:
        name = model.get("sequence_id")
        if name not in wanted:
            continue
        label = model.get("species") or model.get("organism")
        if isinstance(label, str) and label.strip():
            by_id[name].add(label.strip())
        accession = model.get("source_accession")
        if isinstance(accession, str):
            accessions[name].add(re.sub(r"\.\d+$", "", accession.strip()))
    for name, label in supplied.items():
        if name in wanted:
            by_id[name].add(label)
    reviewed = {}
    content = files("pocket_workbench").joinpath(
        "content", "BRCA1_LINKED_LANDMARKS.tsv"
    )
    for line in content.read_text(encoding="utf-8").splitlines()[1:]:
        digest, accession, *_ = line.split("\t")
        if accession in REVIEWED_TAXA:
            reviewed[digest, accession] = REVIEWED_TAXA[accession]
    for name in wanted:
        sequence = data["sequences"].get(name)
        if sequence and not by_id[name]:
            digest = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
            by_id[name].update(
                reviewed[digest, accession]
                for accession in accessions[name]
                if (digest, accession) in reviewed
            )
    return {
        name: next(iter(by_id[name]))
        if len(by_id[name]) == 1
        else CONFLICT
        if by_id[name]
        else UNKNOWN
        for name in sorted(wanted)
    }


def species_rows(*, assignments: dict[str, str], identifiers: list[str]) -> list[dict]:
    """Count each representative ID once, including unknown or conflicting labels."""
    counts = Counter(assignments.get(name, UNKNOWN) for name in set(identifiers))
    return [
        {"species": label, "sequences": count}
        for label, count in sorted(counts.items(), key=lambda x: (-x[1], x[0]))
    ]


def tree_subgroups(
    *, newick: str, identifiers: list[str], count: int
) -> dict[str, str]:
    """Cut the displayed tree into deterministic display clades by terminal count.

    The existing display root is used. The largest splittable clade is split
    first. Multifurcations can make the output contain more groups than asked.
    This is an exploratory tree view, not orthology inference.

    Args:
        newick: Exact saved tree displayed in the app.
        identifiers: Exact saved alignment identifiers.
        count: Target number of display groups.

    Returns:
        Representative ID to display clade label.
    """
    if not 2 <= count <= min(20, len(identifiers)) or len(newick) > 5_000_000:
        raise WorkbenchError("Choose 2–20 display groups within the alignment size")
    try:
        root = Phylo.read(io.StringIO(newick), "newick").root
        tips = [tip.name for tip in root.get_terminals()]
        if set(tips) != set(identifiers) or len(tips) != len(set(tips)):
            raise ValueError("Tree/alignment membership differs")
    except (ValueError, IndexError, TypeError) as error:
        raise WorkbenchError(
            "The saved tree does not match the alignment IDs"
        ) from error
    groups = [root]
    while len(groups) < count:
        candidates = [node for node in groups if node.clades]
        if not candidates:
            break
        target = min(
            candidates,
            key=lambda node: (
                -len(node.get_terminals()),
                min(tip.name for tip in node.get_terminals()),
            ),
        )
        groups.remove(target)
        groups.extend(target.clades)
    return {
        tip.name: f"Clade {i}"
        for i, node in enumerate(
            sorted(
                groups, key=lambda node: min(tip.name for tip in node.get_terminals())
            ),
            1,
        )
        for tip in node.get_terminals()
    }


def subgroup_rows(
    *, assignments: dict[str, str], membership: dict[str, str], identifiers: list[str]
) -> list[dict]:
    """Count unique representative IDs per species and supplied subgroup."""
    counts = Counter(
        (membership.get(name, "Unassigned"), assignments.get(name, UNKNOWN))
        for name in set(identifiers)
    )
    return [
        {"subcluster": group, "species": label, "sequences": count}
        for (group, label), count in sorted(counts.items())
    ]
