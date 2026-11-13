"""Export atom-level pocket selections in the same frame as saved structures."""

from __future__ import annotations

import hashlib
from collections import defaultdict
from pathlib import Path

import gemmi

from .core import WorkbenchError
from .external_tools import (
    EMERALD_URL,
    MOLSTAR_URL,
    PROTEIN_IMAGER_URL,
    RESIDUE_COLUMNS,
    model_file,
    pocket_tables,
    table_tsv,
)


def imager_files(
    *, root: Path, data: dict, selected: list[dict], pocket_ids: dict[str, int]
) -> dict[str, bytes]:
    """Build full protein and selected-lining mmCIFs for one model or one saved fit.

    Args:
        root: Finished project directory.
        data: Original model and pocket records.
        selected: Models in export order, with the saved aligned path for a mobile.
        pocket_ids: Native fpocket IDs keyed by the exact model IDs.

    Returns:
        ZIP member names and bytes, including coverage, hashes and upload guidance.
    """
    ids = [model["sample_id"] for model in selected]
    if len(ids) not in {1, 2} or set(ids) != pocket_ids.keys():
        raise WorkbenchError("Choose one exact fpocket cavity for each exported model")
    table = pocket_tables(root=root, data=data, model_ids=ids)
    wanted = {
        (row["sample_id"], row["pocket_id"])
        for row in table["pockets"]
        if row["method"] == "fpocket"
        and pocket_ids.get(row["sample_id"]) == row["pocket_id"]
    }
    if wanted != set(pocket_ids.items()):
        raise WorkbenchError(
            "A selected native fpocket cavity is absent from the resource"
        )
    membership = [
        row
        for row in table["residues"]
        if row["method"] == "fpocket" and (row["sample_id"], row["pocket_id"]) in wanted
    ]
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in membership:
        grouped[row["sample_id"]].append(row)
    content = {}
    audit = []
    roles = ("reference", "aligned_mobile") if len(ids) == 2 else ("model",)
    for model, role in zip(selected, roles, strict=True):
        sample = model["sample_id"]
        pocket = pocket_ids[sample]
        path = f"protein_imager/{role}_protein.cif"
        content[path] = model_file(root=root, model=model, file_format="cif")
        structure = gemmi.make_structure_from_block(
            gemmi.cif.read_string(content[path].decode("utf-8")).sole_block()
        )
        keys = {
            (row["tool_chain"], row["tool_residue_number"])
            for row in grouped[sample]
            if row["mapping_status"] == "mapped"
        }
        found = set()
        for chain in structure[0]:
            for index in range(len(chain) - 1, -1, -1):
                residue = chain[index]
                key = (chain.name, residue.seqid.num)
                if key in keys and not residue.seqid.icode.strip():
                    found.add(key)
                else:
                    del chain[index]
        pocket_file = ""
        if found:
            structure.remove_empty_chains()
            pocket_file = f"protein_imager/{role}_fpocket_{pocket}.cif"
            content[pocket_file] = (
                structure.make_mmcif_document().as_string().encode("utf-8")
            )
        audit.append(
            {
                "role": role,
                "sample_id": sample,
                "method": "fpocket",
                "pocket_id": pocket,
                "protein_file": path,
                "pocket_file": pocket_file,
                "n_recorded_lining_residues": len(grouped[sample]),
                "n_exported_lining_residues": len(found),
                "n_unexported_lining_residues": len(grouped[sample]) - len(found),
                "protein_sha256": hashlib.sha256(content[path]).hexdigest(),
                "pocket_sha256": hashlib.sha256(content[pocket_file]).hexdigest()
                if pocket_file
                else "",
            }
        )
    content["protein_imager/files.tsv"] = table_tsv(rows=audit, columns=list(audit[0]))
    content["protein_imager/selected_pocket_residues.tsv"] = table_tsv(
        rows=membership, columns=RESIDUE_COLUMNS
    )
    content["protein_imager/START_HERE.txt"] = (
        "Load the supplied proteins and their matching pocket files together.\n"
        "For a pair, the mobile is already superposed on the reference.\n"
        "Keep the saved frame. Each pocket file contains the original atoms of its\n"
        "mapped lining residues. Choose styles and colours in the viewer.\n"
        "files.tsv lists exact IDs, coverage and checksums. A blank pocket path\n"
        "means there were no mapped lining atoms.\n"
        f"\nProtein Imager: {PROTEIN_IMAGER_URL}\n"
        f"Mol* local-file viewer: {MOLSTAR_URL}\n"
        f"EMERALD (load model_inputs.fasta manually): {EMERALD_URL}\n"
    ).encode("utf-8")
    return content
