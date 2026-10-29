# Protein Pocket Workbench 1.4.0

Browse predicted protein pockets, saved structural alignments, pocket-marked
FASTA alignments and trees. The analysis runs in batch on HPC; colleagues can
open the finished resource on their own computer. No cluster account or pocket
predictor is needed for viewing.

## Updating the viewer you already use



Extract the software ZIP into Downloads, unzip then:

```bash
cd "$HOME/Downloads/protein_workbench"

conda env create --file environment.viewer.yml
conda activate pocket_workbench


python -m pip install --upgrade --editable .
pocket-workbench --version

```

The version should be `1.4.0`. Keep this source directory in place because the
installation is editable. If you do not already have the environment, replace
the activation step above with these two commands:

```bash
conda env create --file environment.viewer.yml
conda activate pocket_workbench
```

## Open your existing BRCA1 results

On the Mac used for viewing:

```bash
MADDY_LOCAL="PATH_TO_RESULTS"
pocket-workbench app \
  --resource "$MADDY_LOCAL/results/brca1_workbench_v1_1" \
  --port 8510
```

Open <http://127.0.0.1:8510>. Adjust `MADDY_LOCAL` to your actual local folder.
For the completed quality rerun, use its finished `04_workbench` resource
instead. Resume outputs may be numbered: consult `latest_outputs.json` for
the exact path and wait until the run has finished before copying it.

For any other project:

```bash
pocket-workbench app --resource /path/to/finished_resource --port 8510
```

Use the directory containing **project.json**, prepared models and
**pocket_workbench.duckdb**. Keep the complete resource together. An ordinary
HTML report or DuckDB file alone does not supply the structures needed for 3D.
A review ZIP contains a nested `app_resource.zip`; extract that too for the app.
The app starts in browsing mode, with analysis submission controls hidden.

## Share the updated offline viewer

The app immediately uses the new viewer when opening an existing native
resource. Previously exported HTML files retain the viewer embedded when they
were created. To regenerate a finished offline review from saved evidence:

```bash
pocket-workbench review \
  --resource /path/to/finished_resource \
  --output /path/to/new_review_v1_4
```

Choose a new output folder outside the resource. This exports existing evidence;
it does not run predictors, structural fitting, FAMSA or tree inference. It
creates the review folder and a sibling ZIP, including HTML, formatted Excel,
TSV, DuckDB and `app_resource.zip`. All-pairs reviews can be large: allow space
for the folder and ZIP. Share the ZIP; reviewers extract the whole folder and
open `index.html`. External source/Inspector links require internet.

## Batch analysis and methods

The quality-preparation/prediction workflow and Slurm scripts are retained.
Use [QUALITY_RERUN.md](docs/QUALITY_RERUN.md) for a new quality run and
[HPC_WORKFLOW.md](docs/HPC_WORKFLOW.md) to build from existing predictions.
Do not replace the running 1.2.0 workflow as part of this viewer update.

Keep original full-length predictions as the baseline. Automatic confidence
regions are candidates, not confirmed biological domains; cropping changes
surface geometry and can create exposed sites. Interpret pockets with local
confidence, cut-boundary flags, fit coverage and biological context.

[Methods](docs/METHODS.md) · [Glossary](docs/GLOSSARY.tsv) ·
[Deployment](docs/DEPLOYMENT.md) · [Inputs](docs/INPUTS.md) ·
[DuckDB schema](docs/SCHEMA.md) · [Validation](docs/VALIDATION.md)

Development validation uses explicitly synthetic coordinates and saved evidence.
The real BRCA1 coordinate resource and running Dundee jobs were not available
in this workspace; no new biological predictions or cluster submissions were made.
