# KaroSpaceR

`KaroSpaceR` is a separate, repo-local implementation for exporting KaroSpace-like viewers from `.rds` inputs without modifying the `KaroSpace` repository.

Current scope:

- Pure-R payload normalization and HTML export
- Vendored viewer shell snapshot stored in this repo
- Render parity first: same viewer shell and core interactions
- Marker-gene export for categorical colors
- Optional contact-conditioned interaction markers on top of spatial neighbors
- Automatic sparse gene encoding and packed section arrays for large exports
- Optional `lightweight` export mode for smaller HTML output
- Supported inputs today:
  - `Seurat`
  - Plain R lists with `obs`, `coordinates`, and optional `expression` / `umap`
  - `SingleCellExperiment`
  - `SpatialExperiment`

## Repo workflow

1. Refresh the vendored viewer shell from the reference KaroSpace repo:

```bash
python3 scripts/sync_viewer_shell.py
```

2. Export a viewer from an `.rds` input:

```bash
Rscript scripts/karospace_build_r.R \
  --input input.rds \
  --output viewer.html \
  --groupby sample_id \
  --initial-color cell_type \
  --additional-colors course \
  --assay SCT \
  --top-genes 200 \
  --marker-genes-groupby auto
```

## Config-File Workflow

The recommended backend contract for future automation or a desktop wrapper is a JSON build config.

The recommended discovery contract is a machine-readable inspect command:

```bash
Rscript scripts/karospace_inspect_r.R \
  --input path/to/object.rds \
  --output inspect.json
```

That JSON includes:

- detected defaults such as `groupby`, `initial_color`, `additional_colors`, and `assay`
- available assays and Seurat graph names
- categorical and numeric metadata columns
- gene counts and a preview
- a `default_config` object that can be written directly to the build step

For a future desktop app, the intended flow is:

1. call `karospace_inspect_r.R`
2. let the user edit the returned `default_config`
3. call `karospace_build_r.R --config build.json`

Generate a resolved config from a real dataset with auto-detected defaults:

```bash
Rscript scripts/example_export.R \
  --input path/to/object.rds \
  --write-config build.json \
  --inspect
```

That writes a versioned JSON config with explicit export settings.

Run the export from that config:

```bash
Rscript scripts/karospace_build_r.R --config build.json
```

You can still override individual fields at run time:

```bash
Rscript scripts/karospace_build_r.R \
  --config build.json \
  --output viewer_alt.html \
  --top-genes 300 \
  --write-config build_resolved.json
```

Config fields use snake_case and map directly to the exporter API, for example:

```json
{
  "version": 1,
  "input": "/path/to/object.rds",
  "output": "/path/to/viewer.html",
  "groupby": "orig.ident",
  "initial_color": "seurat_clusters",
  "additional_colors": ["integrated_snn_res.0.47"],
  "assay": "SCT",
  "top_genes_n": 200,
  "neighbor_mode": "spatial",
  "marker_genes_groupby": ["auto"],
  "interaction_markers_groupby": ["seurat_clusters"],
  "title": "IntegratedHeartST"
}
```

## Fastest Real-Dataset Test

If you just want to try one of your real `.rds` files with minimal setup, run:

```bash
Rscript scripts/example_export.R --input path/to/object.rds
```

That script will:

- inspect the object
- auto-pick a likely `groupby`
- auto-pick a likely initial color column
- choose a few extra color columns when possible
- write `<input_stem>_karospacer.html` beside the `.rds`

If you only want to inspect what it would choose:

```bash
Rscript scripts/example_export.R --input path/to/object.rds --inspect
```

If you want to inspect available genes for the selected assay before exporting:

```bash
Rscript scripts/example_export.R \
  --input path/to/object.rds \
  --assay SCT \
  --inspect-genes
```

To filter the gene list:

```bash
Rscript scripts/example_export.R \
  --input path/to/object.rds \
  --assay SCT \
  --inspect-genes \
  --gene-query COL \
  --gene-limit 25
```

If auto-detection picks the wrong columns, override them explicitly:

```bash
Rscript scripts/example_export.R \
  --input path/to/object.rds \
  --groupby sample_id \
  --initial-color cell_type \
  --additional-colors course,condition \
  --assay SCT \
  --genes CXCL8,COL1A1 \
  --output viewer.html
```

If you want the exporter to preload the highest-expressed genes from the chosen assay instead of naming genes manually:

```bash
Rscript scripts/example_export.R \
  --input path/to/object.rds \
  --assay SCT \
  --top-genes 200
```

`--top-genes 200` selects the top 200 genes by mean expression across all exported cells/spots in the selected assay.

If you want a smaller standalone HTML and do not need precomputed marker or interaction panels, use lightweight mode:

```bash
Rscript scripts/example_export.R \
  --input path/to/object.rds \
  --assay SCT \
  --top-genes 200 \
  --lightweight
```

Phase 5 behavior:

- large section arrays are packed automatically into base64 float32/uint32 payloads
- mostly-zero gene vectors are exported as sparse index/value arrays
- `--lightweight` disables heavy precomputed analyses unless you explicitly turn them back on

To precompute contact-conditioned interaction markers for a categorical color:

```bash
Rscript scripts/example_export.R \
  --input path/to/object.rds \
  --groupby sample_id \
  --initial-color cell_type \
  --additional-colors course \
  --neighbor-mode spatial \
  --marker-genes-groupby auto \
  --interaction-markers-groupby cell_type \
  --interaction-markers-top-targets 8 \
  --interaction-markers-top-genes 12 \
  --interaction-markers-min-cells 30 \
  --interaction-markers-min-neighbors 1 \
  --output viewer.html
```

Marker genes default to the exported categorical colors in `example_export.R`. Interaction markers stay opt-in because they are more expensive and depend on the neighbor graph.

## Plain-list input contract

The simplest supported `.rds` object is a list with:

- `obs`: data frame with one row per cell
- `coordinates`: numeric matrix/data frame with two columns
- `expression`: optional genes x cells matrix-like object
- `umap`: optional numeric matrix/data frame with two columns
- `gene_names`: optional character vector matching the expression rows
