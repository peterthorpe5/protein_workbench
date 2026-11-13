# Motif search, 3D snapshots and reference inputs

## Find an amino acid string in the model

Open **Models & pockets**, select one model, and use **Amino acid search**.
For the intact dog `NP_001013434` model 0, select native fpocket **33** and search
for `ESQDR`. The hit selector shows **1757–1761** and **1762–1766** on the
verified representative axis. Choose either range or **All occurrences**.
The colour picker changes the motif highlight and its Nightingale track.
Pocket lining, the selected BRCT window and the individual residue selection
remain separate. Where the motif and pocket overlap, the 3D view marks the
overlap while keeping the pocket colour visible.

Search is literal and case insensitive, allows 1–100 amino acid letters,
ignores whitespace and includes overlapping hits. `X` matches only `X`.
The viewer prefers the verified representative sequence; otherwise it labels
the axis as prepared-model numbering. A trimmed model can have an occurrence
on the representative sequence with no coordinates. Only mapped residues are
highlighted in 3D.

In **Compare structures**, search the reference, mobile or both models in the
saved frame. Recorded landmarks are available for either model when their
sequence matches a reviewed sequence. The two model axes are checked separately.
No new fit or pocket prediction is run. Downloaded 3D HTML has the same tools;
an older HTML export needs reopening in the updated app and exporting again.

## Download the view

Rotate or zoom to the view you want, then select **Download PDF** on the 3D
toolbar. This saves a static A4 landscape page with the rendered scene, model
IDs and a legend for selected pockets, landmarks, motif hits and picked residue.
It works in the offline HTML and preserves the active colours. The PNG button
still provides the raw scene image.

## Missing structures and sequences

Use **Download a missing model or sequence** on **Models & pockets**, for a
missing member of a saved alignment, or in **HOG / model groups**. Requests run
only when you click Fetch. Public coordinates and FASTA download separately
with their accession, source, retrieval date and checksum.

- **AlphaFold DB:** supply a canonical UniProt accession. The app checks the
  returned model against its declared sequence and any exact saved full sequence.
  **Preview the public model with Nightingale** opens it separately for motif
  search, residue selection and scene export. Its original mmCIF remains available.
  It does not infer a UniProt accession from a RefSeq ID or choose between
  several fragments. To recover one of Maddy's original predictions, import
  that exact de novo model through **Start / import**.
- **NCBI Protein:** supply a full RefSeq or GenBank protein accession, optionally
  with its version. The returned record must match. The app reports whether
  its sequence agrees with the saved full sequence, when available.
  Import the FASTA as input to a new analysis when you want to use it. The saved
  alignment is unchanged. Hosts can set `POCKET_WORKBENCH_NCBI_EMAIL` for the
  NCBI contact parameter.

These downloads do not attach public coordinates to an existing pocket record.
See [AlphaFold DB API documentation](https://alphafold.ebi.ac.uk/api-doc) and
the [NCBI EFetch formats](https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/).

## Species and subcluster composition

Open **Species composition** on a completed alignment or HOG/model group.
The table, chart and TSV count each representative ID once. Distinct IDs with
identical amino acid strings remain distinct. Species come from explicit
metadata or the reviewed BRCA1 dog/human accession and exact sequence. Add a
tab-separated
`sequence_id` and `species` file for the remaining records. Missing or
conflicting labels have their own categories.

For **Species by subcluster**, supply a `sequence_id` and `subcluster` TSV.
Each ID may have one assignment. Unassigned members stay in the summary.
On a saved alignment, **Preview a cut of the saved tree** can show an
exploratory composition by display clade. It splits the largest remaining
branch by terminal count; a multifurcation can yield more than the target
number of groups. These display groups do not alter the HOG, infer orthology
or change saved analyses.
Choose a subcluster to compare its species counts with the complete selection.
The tables and TSVs retain all species, even when the charts show a subset.
