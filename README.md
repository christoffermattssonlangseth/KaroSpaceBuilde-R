# KaroSpaceBuildR

`KaroSpaceBuildR` is a separate, repo-local implementation for exporting KaroSpace-like viewers from `.rds` inputs without modifying the `KaroSpace` repository.

Current scope:

- Pure-R payload normalization and HTML export
- Vendored viewer shell snapshot stored in this repo
- Render parity first: same viewer shell and core interactions
- Supported inputs today:
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
  --additional-colors course
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
- write `<input_stem>_karospace_buildr.html` beside the `.rds`

If you only want to inspect what it would choose:

```bash
Rscript scripts/example_export.R --input path/to/object.rds --inspect
```

If auto-detection picks the wrong columns, override them explicitly:

```bash
Rscript scripts/example_export.R \
  --input path/to/object.rds \
  --groupby sample_id \
  --initial-color cell_type \
  --additional-colors course,condition \
  --output viewer.html
```

## Plain-list input contract

The simplest supported `.rds` object is a list with:

- `obs`: data frame with one row per cell
- `coordinates`: numeric matrix/data frame with two columns
- `expression`: optional genes x cells matrix-like object
- `umap`: optional numeric matrix/data frame with two columns
- `gene_names`: optional character vector matching the expression rows
