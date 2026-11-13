"""Viewer controls for public reference downloads and species composition."""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from .core import WorkbenchError
from .external_tools import canonical_uniprot, table_tsv
from .public_records import (
    fetch_afdb_model,
    fetch_ncbi_fasta,
    ncbi_protein_accession,
    public_model_document,
)
from .species import (
    CONFLICT,
    UNKNOWN,
    read_species_map,
    species_assignments,
    species_rows,
    subgroup_rows,
    tree_subgroups,
)

LOGGER = logging.getLogger(__name__)


def public_record_panel(
    *,
    key: str,
    root: Path,
    data: dict,
    model: dict | None = None,
    sequence_id: str = "",
    sequences: dict[str, str] | None = None,
) -> None:
    """Offer explicit public downloads, bound to the current selection and project.

    Args:
        key: Unique page-local control prefix.
        root: Open project directory.
        data: Saved resource records, used for exact sequence checks.
        model: Selected local model, when applicable.
        sequence_id: Selected representative ID on the sequence page.
        sequences: Optional exact saved alignment sequences keyed by ID.
    """
    model = model or {}
    sequence_id = sequence_id or model.get("sequence_id") or ""
    expected = (sequences or {}).get(sequence_id) or data["sequences"].get(sequence_id)
    saved = bool(model.get("prepared_model")) or bool(
        sequence_id
        and any(
            row.get("sequence_id") == sequence_id and row.get("prepared_model")
            for row in data["models"]
        )
    )
    exact_model_seq = model.get("sequence")
    if not expected and model.get("input_scope", "full_model") == "full_model":
        expected = exact_model_seq
    identifiers = [model.get("source_accession"), sequence_id, model.get("sample_id")]
    ncbi = next(
        (
            value
            for source in identifiers
            if (value := ncbi_protein_accession(value=source))
        ),
        "",
    )
    uniprot = next(
        (value for source in identifiers if (value := canonical_uniprot(value=source))),
        "",
    )
    context = (str(root), data.get("revision"), model.get("sample_id"), sequence_id)
    tag = hashlib.sha256(repr(context).encode()).hexdigest()[:12]
    with st.expander("Download a missing model or sequence"):
        st.caption(
            "Download public reference files by exact accession. This does not attach "
            "coordinates or sequences to the saved project or carry over its pockets."
        )
        providers = [("NCBI protein FASTA", ncbi, fetch_ncbi_fasta)]
        if not saved:
            providers.insert(0, ("AlphaFold DB mmCIF", uniprot, fetch_afdb_model))
        for label, preset, fetch in providers:
            provider_key = "afdb" if label.startswith("AlphaFold") else "ncbi"
            value = (
                st.text_input(
                    "UniProt accession for AlphaFold DB"
                    if provider_key == "afdb"
                    else "NCBI protein accession",
                    value=preset,
                    key=f"{key}_{tag}_{provider_key}_id",
                    max_chars=30,
                    help="Use the full accession. NCBI versions are optional.",
                )
                .strip()
                .upper()
            )
            token = (context, value)
            state_key = f"{key}_{provider_key}_download"
            state = st.session_state.get(state_key)
            if state and state["token"] != token:
                st.session_state.pop(state_key, None)
            if st.button(
                f"Fetch {label}", key=f"{key}_{provider_key}_fetch", disabled=not value
            ):
                st.session_state.pop(state_key, None)
                try:
                    opts = {"accession": value}
                    if provider_key == "afdb":
                        opts["expected"] = expected
                    with st.spinner(f"Fetching {label}…"):
                        record, payload = fetch(**opts)
                    st.session_state[state_key] = dict(
                        token=token, record=record, payload=payload
                    )
                except (WorkbenchError, OSError) as error:
                    LOGGER.warning(
                        "Reference lookup failed for %s: %s", provider_key, error
                    )
                    st.warning(str(error))
            if state := st.session_state.get(state_key):
                record = state["record"]
                accession = record["accession"]
                suffix = "cif" if provider_key == "afdb" else "fasta"
                if expected:
                    identity = (
                        "Exact saved sequence match"
                        if expected == record["sequence"]
                        else "Differs from the saved sequence"
                    )
                else:
                    identity = "No saved full sequence was available to compare"
                species = f" · {record['species']}" if record.get("species") else ""
                st.caption(
                    f"{accession}{species} · {len(record['sequence'])} residues. "
                    f"{identity}."
                )
                st.download_button(
                    f"Download {label}",
                    state["payload"],
                    file_name=f"{accession}.{suffix}",
                    mime="chemical/x-mmcif" if suffix == "cif" else "text/plain",
                    key=f"{key}_{provider_key}_file",
                    on_click="ignore",
                )
                if provider_key == "afdb" and st.checkbox(
                    "Preview the public model with Nightingale",
                    key=f"{key}_{tag}_{accession}_preview",
                    help="Opens a separate view without project pocket predictions.",
                ):
                    try:
                        document = state.get("document")
                        if document is None:
                            document = public_model_document(
                                record=record, payload=state["payload"]
                            )
                            state["document"] = document
                        st.iframe(
                            document,
                            height="content",
                            width="stretch",
                            alt="Public AlphaFold model with linked amino acid search",
                        )
                        st.download_button(
                            "Download offline public model viewer",
                            document,
                            file_name=f"{accession}_public_model.html",
                            mime="text/html",
                            on_click="ignore",
                        )
                    except (WorkbenchError, OSError, RuntimeError) as error:
                        st.warning(
                            f"The model download is ready; preview unavailable: {error}"
                        )
                with st.expander(f"{label} provenance"):
                    audit = {k: v for k, v in record.items() if k != "sequence"}
                    st.json(audit)
                    st.download_button(
                        f"Download {label} provenance (TSV)",
                        table_tsv(rows=[audit], columns=list(audit)),
                        file_name=f"{accession}_{provider_key}_source.tsv",
                        mime="text/tab-separated-values",
                        on_click="ignore",
                    )
    if not saved:
        st.caption("Add a missing original de novo model through Start / import.")


def composition_panel(
    *,
    data: dict,
    identifiers: list[str],
    key: str,
    newick: str = "",
    show_groups: bool = False,
) -> None:
    """Display sequence-level composition for the selected group or alignment.

    Args:
        data: Open workbench project.
        identifiers: Representative sequence IDs in the current selection.
        key: Unique widget and export prefix.
        newick: Matching saved display tree, when available.
        show_groups: Enable supplied subcluster or exploratory tree breakdown.
    """
    with st.expander("Species composition", expanded=True):
        if not identifiers:
            st.info("No representative sequence IDs are available for this selection.")
            return
        st.caption(
            "Each representative sequence ID is counted once, even with several models."
        )
        supplied = {}
        metadata = st.file_uploader(
            "Optional species TSV",
            type=["tsv"],
            key=f"{key}_species_metadata",
            help="Two columns: sequence_id and species. Match the saved FASTA IDs.",
        )
        if metadata:
            try:
                if metadata.size > 1_000_000:
                    raise WorkbenchError("The species TSV exceeds the 1 MB limit")
                supplied = read_species_map(text=metadata.getvalue().decode("utf-8"))
            except (WorkbenchError, UnicodeError) as error:
                st.warning(str(error))
        assignments = species_assignments(
            data=data, identifiers=identifiers, supplied=supplied
        )
        rows = species_rows(assignments=assignments, identifiers=identifiers)
        st.dataframe(rows, hide_index=True, width="stretch")
        if all(row["species"] in {UNKNOWN, CONFLICT} for row in rows):
            st.info(
                "Species are unresolved. Add a TSV with exact sequence IDs and species."
            )
        shown = rows[:25]
        if len(shown) < len(rows):
            st.caption(
                "The chart shows the 25 largest counts. The table and TSV are complete."
            )
        st.plotly_chart(
            px.bar(
                pd.DataFrame(shown[::-1]),
                x="sequences",
                y="species",
                orientation="h",
                height=max(230, 24 * len(shown) + 105),
                labels={"sequences": "Representative sequences", "species": "Species"},
            ).update_xaxes(
                dtick=1 if max(row["sequences"] for row in shown) <= 15 else None
            ),
            key=f"{key}_species_plot",
            width="stretch",
        )
        st.download_button(
            "Download species counts (TSV)",
            table_tsv(rows=rows, columns=("species", "sequences")),
            file_name=f"{key}_species_counts.tsv",
            key=f"{key}_species_export",
        )
        if not show_groups:
            return
        st.markdown("**Species by subcluster**")
        membership = {}
        uploaded = st.file_uploader(
            "Optional subcluster TSV",
            type=["tsv"],
            key=f"{key}_cluster_metadata",
            help="Two columns: sequence_id and subcluster. Match the displayed IDs.",
        )
        if uploaded:
            try:
                if uploaded.size > 1_000_000:
                    raise WorkbenchError("The subcluster TSV exceeds the 1 MB limit")
                membership = read_species_map(
                    text=uploaded.getvalue().decode("utf-8"), label_column="subcluster"
                )
            except (WorkbenchError, UnicodeError) as error:
                st.warning(str(error))
        elif newick and len(set(identifiers)) >= 3:
            if st.checkbox(
                "Preview a cut of the saved tree", key=f"{key}_tree_preview"
            ):
                count = st.slider(
                    "Target display clades",
                    2,
                    min(20, len(set(identifiers))),
                    2,
                    key=f"{key}_tree_count",
                )
                try:
                    membership = tree_subgroups(
                        newick=newick, identifiers=identifiers, count=count
                    )
                    st.caption(
                        "Exploratory split using the saved tree's display root. "
                        "It does not change group membership or infer orthology."
                    )
                except WorkbenchError as error:
                    st.warning(str(error))
        if not membership:
            return
        clusters = subgroup_rows(
            assignments=assignments, membership=membership, identifiers=identifiers
        )
        st.dataframe(clusters, hide_index=True, width="stretch")
        options = sorted({row["subcluster"] for row in clusters})
        selected = st.selectbox(
            "Subcluster to compare with the complete selection",
            options,
            key=f"{key}_subcluster_compare",
        )
        totals = {row["species"]: row["sequences"] for row in rows}
        counts = {
            row["species"]: row["sequences"]
            for row in clusters
            if row["subcluster"] == selected
        }
        labels = sorted(
            totals, key=lambda label: (-counts.get(label, 0), -totals[label], label)
        )[:25]
        chart = pd.DataFrame(
            [
                {
                    "species": label,
                    "Complete selection": totals[label],
                    f"Subcluster: {selected}": counts.get(label, 0),
                }
                for label in labels[::-1]
            ]
        )
        st.plotly_chart(
            px.bar(
                chart,
                y="species",
                x=["Complete selection", f"Subcluster: {selected}"],
                orientation="h",
                barmode="group",
                height=max(250, 32 * len(labels) + 95),
                labels={
                    "value": "Representative sequences",
                    "variable": "Selection",
                    "species": "Species",
                },
            ),
            key=f"{key}_subcluster_plot",
            width="stretch",
        )
        st.download_button(
            "Download subcluster species counts (TSV)",
            table_tsv(rows=clusters, columns=("subcluster", "species", "sequences")),
            file_name=f"{key}_subcluster_species.tsv",
            key=f"{key}_subcluster_export",
        )
