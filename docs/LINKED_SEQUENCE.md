# Linked structure and pockets

Use **Models & pockets → Linked structure and pockets** on a completed native
resource with prepared coordinates. Select one model, then choose **fpocket**,
**P2Rank**, or its rescoring evidence where available. You can do this with an
existing resource; no new calculations or DuckDB rebuild are needed.

1. **Find a pocket.** Type a native pocket number such as `33`, or search by
   `rank 3`, in the filter above the lining tracks. Click a track or its label to
   select that pocket in the existing 3D view and focus the sequence window.
   Short segments represent lining residues; gaps remain visible. A pocket ID
   identifies a particular prediction. The review rank is its order within the
   model and is shown separately. Track labels report how many lining residues
   map to the displayed axis.
2. **Find a residue.** Click a sequence letter or the confidence trace, or type
   its position and press **Find residue**. The selected residue is marked in
   teal in 3D. Picking an atom in the model highlights the corresponding track
   position. The label gives its original chain and residue number, with an
   insertion code where present, and pLDDT if recorded. Clear the residue and
   landmark highlights without changing the pocket using **Clear highlights**.
3. **Navigate and share.** Drag the overview window, select **Full sequence**
   or return to the selected pocket with **Show selected pocket region**. Use
   **Download offline model pocket page** to share the interactive HTML for the
   selected model and evidence. The page works locally without an app install or
   network request. A browser with WebGL is needed for 3D.

### Maddy's BRCA1 example

Choose the **intact** `fold_np_001013434_03112026_model_0` and **fpocket**. If the
model is missing, open the completed untrimmed resource that includes the dog
models. Search for pocket **33**, then click **ESQDR copy 2** or enter **1766**
in the residue box. The motif and the pocket can stay highlighted together.
The motif uses magenta, the selected residue uses teal, and the pocket keeps
its existing rank colour. The BRCT1 and BRCT2 windows mark the regions used for
the comparison. Model 1 pocket 9, model 3 pocket 41 and model 4 pocket 34 are
the other candidates identified in the completed review. Use the evidence
pack for how those candidates were compared.

Domain windows are violet by default and motifs magenta. To change them, open
**Domain and motif colours** alongside the model. Under **Pocket appearance**,
choose **Choose pocket colours** for a custom lining colour. The default stays
rank-based for one model and blue/orange for a comparison. A comparison can
use a different custom pocket colour for each model. These controls affect
both the lining display and the matching track colour where shown.

### Coordinate and evidence scope

The viewer checks that the imported residue map agrees with the selected
representative FASTA. An exactly mapped cropped region uses the parent FASTA
axis, but a position outside the crop reports that there are no coordinates.
If an exact FASTA mapping cannot be checked, the track explicitly uses
prepared-model positions. The displayed original residue identifier is shown
separately. A model without a usable residue map retains the existing 3D view.

The BRCT and ESQDR landmarks come from the three representative sequences in
the BRCT review. They are drawn only when the displayed sequence has the exact
reviewed SHA-256 fingerprint. They are comparison windows and sequence motifs;
the viewer does not infer domains in other sequences. Changing the display
does not alter the predictions. A highlight represents a predicted pocket's
lining, rather than a validated binding site.

**Compare structures** retains its saved two-model 3D overlay and comparison
metrics. The linked sequence panel uses one model at a time, so sequence
coordinates from distinct models are not conflated. See
[external tools](EXTERNAL_TOOLS.md#export-the-current-3d-comparison) to
download that saved overlay and its selected atomic linings.
