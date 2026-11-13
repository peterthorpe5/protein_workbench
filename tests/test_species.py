"""Test explicit species counts, missing labels and exploratory group membership."""

from hashlib import sha256

import pytest

from pocket_workbench import species
from pocket_workbench.core import WorkbenchError
from pocket_workbench.species import (
    CONFLICT,
    UNKNOWN,
    read_species_map,
    species_assignments,
    species_rows,
    subgroup_rows,
    tree_subgroups,
)


def test_explicit_metadata_counts_sequence_ids_and_preserves_conflicts():
    """Replicates do not inflate counts, and conflicting assignments remain labelled."""
    data = {
        "sequences": {"A": "ACDE", "B": "ACDE", "C": "AAAA", "D": "CCCC"},
        "models": [
            {"sequence_id": "A", "species": "Human"},
            {"sequence_id": "A", "species": "Human"},
            {"sequence_id": "B", "organism": "Dog"},
            {"sequence_id": "C", "species": "Dog"},
        ],
    }
    ids = ["A", "A", "B", "C", "D"]
    supplied = read_species_map(text="sequence_id\tspecies\nC\tHuman\nZ\tCat\n")
    assignments = species_assignments(data=data, identifiers=ids, supplied=supplied)
    assert assignments == {"A": "Human", "B": "Dog", "C": CONFLICT, "D": UNKNOWN}
    assert (
        sum(
            row["sequences"]
            for row in species_rows(assignments=assignments, identifiers=ids)
        )
        == 4
    )
    membership = read_species_map(
        text="sequence_id\tsubcluster\nA\t01\nB\t01\nC\t02\n", label_column="subcluster"
    )
    rows = subgroup_rows(
        assignments=assignments, membership=membership, identifiers=ids
    )
    assert {row["subcluster"] for row in rows} == {"01", "02", "Unassigned"}
    assert sum(row["sequences"] for row in rows) == 4


def test_reviewed_species_requires_accession_as_well_as_sequence(tmp_path, monkeypatch):
    """An identical sequence from an unrelated record cannot assign its species."""
    landmarks = tmp_path / "content" / "BRCA1_LINKED_LANDMARKS.tsv"
    landmarks.parent.mkdir()
    digest = sha256(b"ACDE").hexdigest()
    landmarks.write_text(
        f"digest\taccession\n{digest}\tNP_001013434\n", encoding="utf-8"
    )
    monkeypatch.setattr(species, "files", lambda _package: tmp_path)
    sequences = {
        "NP_001013434.1": "ACDE",
        "alias": "ACDE",
        "same_sequence": "ACDE",
        "wrong_sequence": "ACDF",
    }
    data = {
        "sequences": sequences,
        "models": [
            {"sequence_id": name, "source_accession": "NP_001013434.1"}
            for name in ("alias", "wrong_sequence")
        ],
    }
    dog = species.REVIEWED_TAXA["NP_001013434"]
    assert species_assignments(data=data, identifiers=list(sequences), supplied={}) == {
        "NP_001013434.1": dog,
        "alias": dog,
        "same_sequence": UNKNOWN,
        "wrong_sequence": UNKNOWN,
    }


@pytest.mark.parametrize(
    "text",
    [
        "",
        "sequence_id,species\nA,Dog",
        "sequence_id\tspecies\nA\t",
        "sequence_id\tspecies\nA\tDog\nA\tDog",
        "sequence_id\tspecies\n",
    ],
)
def test_species_map_refuses_malformed_or_duplicate_ids(text):
    """An ambiguous sidecar is never silently merged."""
    with pytest.raises(WorkbenchError):
        read_species_map(text=text)


def test_exploratory_tree_subgroups_keep_exact_unique_tips():
    """Cuts use the saved display root and report all members of a multifurcation."""
    ids = list("ABCD")
    tree = "((A:1,B:1):2,C:3,D:4);"
    groups = tree_subgroups(newick=tree, identifiers=ids, count=2)
    assert len(set(groups.values())) == 3 and groups["A"] == groups["B"]
    four = tree_subgroups(newick=tree, identifiers=ids, count=4)
    assert set(four) == set(ids) and len(set(four.values())) == 4
    with pytest.raises(WorkbenchError, match="does not match"):
        tree_subgroups(newick=tree, identifiers=list("ABCE"), count=2)
